import sys
import os
import argparse
import shutil
import datetime
import numpy as np
from netCDF4 import Dataset

class PackNetCDF(object):
    def __init__(self):
        self.fin, self.fout, self.overwrite = self.parse_cmd()
    # number of values 2 and 4 byte unsigned integers
    outResolutionShort = 2.0**16 - 2
    outResolutionLong = 2.0**32 - 2  # for unknown reason 2**32 produces wrong results
                                     # try it anyways - hvw

    # coordinate variables to prevent from compression even if they are 2D
    exclude = ('lon', 'lat', 'slon', 'slat', 'slonu', 'slatu', 'slonv',
               'slatv', 'time', 'time_bnds', 'rlon', 'rlat', 'level_bnds',
               'level', 'levels')
    
    def parse_cmd(self):
        parser = argparse.ArgumentParser(description='Compress a netcdf file by ' +
        'first converting float32 and float64 type2D ,3D and 4D fields into ' +
        'integers with offset and scaling factor and then compressing with zlib ' +
        'compression.')
        parser.add_argument('-o', dest='fout',
                            help='compressed netcdf file', metavar='OUTFILE')
        parser.add_argument('-W', default=False, action='store_true',
                            dest='overwrite', help='replace input file, ' +
                            'overrides -o option.')
        parser.add_argument('fin', help='input file', metavar='INFILE')
        args = vars(parser.parse_args())

        # check input file
        fin = os.path.realpath(args['fin'])
        if not os.path.exists(fin):
            parser.error('input file {} does not exist.'.format(fin))
        dir_in = os.path.dirname(fin)

        # check output file
        if args['overwrite']:
            fout = fin + '.tmp'
        else:
            try:
                # output file specified
                fout = os.path.realpath(args['fout'])
                if not os.path.exists(os.path.dirname(fout)):
                    parser.error('path to output file {} does not exist.'
                                 .format(fout))
            except AttributeError:
                # no output file specified
                dir_out = os.path.join(dir_in, 'compress')
                if not os.path.exists(dir_out):
                    print('creating {}'.format(dir_out))
                    os.mkdir(dir_out)
                fout = os.path.join(dir_out, os.path.basename(fin))
        return((fin, fout, args['overwrite']))

        def cp_input(self):
            shutil.copy(self.fin, self.fout)
        
        # def mk_attributes(self):
        #     dsout = Dataset(fout, 'w', format='NETCDF4')
        #     datetime.datetime.now().ctime() + ': ' + 
            
            
            


 

    ############################################################################
    # Open input and output files, then first copy dimensions and global
    # attributes and then copy variables in compressed form if possible
    ############################################################################
bla = '''
    dsin=Dataset(fin, 'r', format='NETCDF4')
    dsout=Dataset(fout, 'w', format='NETCDF4')

    # Copy global attributes and modify or create history attribute
    cdate = datetime.datetime.now()
    histattr = cdate.ctime()+': '+' '.join(sys.argv)

    ehist = False

    for k in dsin.ncattrs():
        attr = dsin.getncattr(k)
        if k == 'history':
            ehist = True
            attr = histattr + '\n' + attr
            # The hstory attribute can not be stored properly as 2D character array with
            # netcdf-python. The following two lines do not help, unfortunately.
            #attr = attr.encode('ascii','ignore')
            #attr = attr.split('\n')
        if type(attr) is unicode:
            dsout.setncattr(k,attr.encode('ascii','ignore'))
        else:
            dsout.setncattr(k,attr)

    if not ehist:
        dsout.setncattr('history',histattr.encode('ascii','ignore'))

    #Copy dimensions
    for dname, the_dim in dsin.dimensions.iteritems():
        dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    # Copy variables
    for v_name, varin in dsin.variables.iteritems():

        vd = varin.dimensions

        if len(vd) >= 2 and not v_name in exclude and \
           (varin.datatype == 'float32' or varin.datatype == 'float64' or \
            varin.datatype == 'uint16' or varin.datatype == 'uint32'):

            # check range, computed offset and scaling, and check if variable is
            # well behaved (short integer ok) or highly skewed (long integer necessary)
            minVal = np.min(varin[:])
            maxVal = np.max(varin[:])
            meanVal = np.mean(varin[:])

            if (meanVal-minVal) >= (maxVal-minVal)/1000.:
                outResolution=outResolutionShort
                intType = 'u2'
            else:
                outResolution=outResolutionLong
                intType = 'u4'

            #outResolution=outResolutionShort
            #intType = 'u2'

            # compress variable
            print 'variable',v_name,varin.datatype,' compressed as ',intType,minVal,maxVal,meanVal
            if len(vd) == 2:
                chunksizes= (len(dsin.dimensions[vd[0]]),len(dsin.dimensions[vd[1]]))
            elif len(vd) == 3:
                chunksizes= (1,len(dsin.dimensions[vd[1]]),len(dsin.dimensions[vd[2]]))
            elif len(vd) == 4:
                chunksizes= (1,1,len(dsin.dimensions[vd[2]]),len(dsin.dimensions[vd[3]]))
            else:
                print('Warning: not tested for 5-D or higher dimension arrays')
                chunksizes= (1,1,1,len(dsin.dimensions[vd[3]]),len(dsin.dimensions[vd[4]]))
            outVar = dsout.createVariable(v_name,intType,varin.dimensions,zlib=True,complevel=9,\
                                          chunksizes=chunksizes)
            add_offset=minVal
            if maxVal == minVal:
                scale_factor = 1.
            else:
                scale_factor=(maxVal-minVal)/outResolution
            outVar.setncattr('scale_factor',scale_factor)
            outVar.setncattr('add_offset',add_offset)
            outVar.set_auto_maskandscale(True)

        else:

            outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
            print 'preserved variable',v_name,varin.datatype

        # Copy remaining variable attributes except fill value
        for k in varin.ncattrs():
            if k != '_FillValue':
                # encode attribute as simple ASCII instead of unicode
                attr = varin.getncattr(k)
                if type(attr) is unicode:
                    outVar.setncattr(k,attr.encode('ascii','ignore'))
                else:
                    outVar.setncattr(k,attr)

        outVar[:] = varin[:]

    # close the input and output files
    dsout.close() 
    dsin.close()


    if overwrite:
        os.rename(fout,fin)
'''

if __name__ == "__main__":
    P = PackNetCDF()
    print (P.fin, P.fout, P.overwrite)
    print(P.outResolutionShort)
