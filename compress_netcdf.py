#!/usr/bin/env python

import sys
import os
import getopt
import datetime
import numpy as np
from netCDF4 import Dataset

####################
# Argument parsing #
####################
 

    fin = ''
    fout = ''
    overwrite = False

    try:
        myopts, args = getopt.getopt(sys.argv[1:],"i:o:W")
    except getopt.GetoptError as e:
        print (str(e))
        print("Usage: %s [-W] -i infile -o outfile" % sys.argv[0])
        sys.exit(2)

    for o, a in myopts:
        if o == '-i':
            fin=a
        elif o == '-o':
            fout=a
        elif o == '-W':
            overwrite = True
        else:
            print("Usage: %s [-W] -i infile -o outfile" % sys.argv[0])
            sys.exit(3)

    if fin == '':
        print("Usage: %s [-W] -i infile -o outfile" % sys.argv[0])
        sys.exit(4)

    # full path for input file
    dir_in=os.path.dirname(fin)
    fname=os.path.basename(fin)
    if dir_in == '':
        dir_in = os.getcwd()
        fin = dir_in+os.sep+fname

    # full path for output file
    if overwrite:
        fout = fin+'.tmp'
    elif fout == '':
        fout=dir_in+os.sep+'compress'+os.sep+fname

    dir_out=os.path.dirname(fout)
    if dir_out == '':
        dir_out = os.getcwd()
        fout = dir_out+os.sep+os.path.basename(fout)

    print ("Input file : %s ; output file: %s" % (fin,fout) )

    # check if input file and output directory exist
    if not os.path.isfile(fin):
        print('input file not found')
        sys.exit(5)

    if not os.path.isdir(dir_out):
        print('generating directory '+dir_out)
        os.mkdir(dir_out)

    #############
    # Constants #
    #############

    # number of values 2 and 4 byte unsigned integers
    outResolutionShort = 2.**16-1.
    outResolutionLong = 2.**31-1. # for unknown reason 2**32 produces wrong results

    # coordinate variables to prevent from compression even if they are 2D
    exclude = ('lon','lat','slon','slat','slonu','slatu','slonv','slatv',\
               'time','time_bnds','rlon','rlat','level_bnds','level','levels')

    ############################################################################
    # Open input and output files, then first copy dimensions and global
    # attributes and then copy variables in compressed form if possible
    ############################################################################

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

if __name__ == "__main__":
    print "running"
