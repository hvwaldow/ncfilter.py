import sys
import os
import argparse
import datetime
from collections import OrderedDict
import numpy as np
from netCDF4 import Dataset


class NcFilter(object):

    def __init__(self, origin):
        '''
        Read all the meta-data of the source file
        into reasonable data-structures.
        '''
        self.origin = origin
        with Dataset(origin, 'r') as dsin:
            # Global attributes
            self.glob_atts = OrderedDict([(x, dsin.getncattr(x))
                                          for x in dsin.ncattrs()])
            # Dimensions
            dim_sizes = [None if x.isunlimited() else len(x)
                         for x in dsin.dimensions.values()]
            self.dims = OrderedDict(zip(dsin.dimensions.keys(), dim_sizes))
            # variables
            # All keys have to be present! In case of no attributes use empty dict.
            self.variables = OrderedDict([(x.name, {
                'dtype': x.dtype,
                'dimensions': x.dimensions,  # tuple
                'attributes': self._get_var_attrs(x)})
                for x in dsin.variables.values()])
        self.newdata = {}

    def _get_var_attrs(self, v):
        return(OrderedDict([(x, v.getncattr(x)) for x in v.ncattrs()]))

    def _get_origin_values(self, varname):
        with Dataset(self.origin, 'r') as ds:
            return(ds.variables[varname][:])

    def _mk_empty_data(self, varname, dimensions, dtyp):
        return({varname: np.ma.MaskedArray(
            np.zeros(dimensions, dtype=np.dtype(dtyp)),
            mask=True)})

    def _get_dimshape(self, varnam):
        dimshape = tuple([self.dims[dimname]
                          for dimname in self.variables[varnam]['dimensions']])
        return(dimshape)

    def write(self, outfile):
        '''
        Creates <outfile> with meta-data as in
          self.glob_atts
          self.dims
          self.variables
        Data in self.newdata ( {<varname>: np.array(...), ...} ) replaces
        the respective data of <varname> in the original file.
        New <varname>s in <self.newdata> have to be present in <self.variables>.
        '''
        dsout = Dataset(outfile, "w")
        dsout.set_auto_mask(True)

        # sanity checks
        if not (type(self.newdata) == dict):
            sys.exit("<self.newdata> has to be a dictionary")
        if not set(self.newdata.keys()) <= set(self.variables.keys()):
            sys.exit("<self.newdata> has undefined variable names: {}"
                     .format(set(self.newdata.keys())
                             .intersection(set(self.variables.keys()))))

        # write global attributes
        dsout.setncatts(self.glob_atts)

        # write dimensions
        for dnam, dsiz in self.dims.iteritems():
            dsout.createDimension(dnam, dsiz)

        # define variables (meta only)
        for vn, v in self.variables.iteritems():
            vout = dsout.createVariable(vn, v['dtype'],
                                        dimensions=v['dimensions'])
            vout.setncatts(v['attributes'])

        # variables to be identically copied (data):
        vcp = set(self.variables.keys()) - set(self.newdata.keys())
        with Dataset(self.origin, "r") as dsin:
            for v in vcp:
                dsout.variables[v][:] = dsin.variables[v][:]

        # variables with new data
        for v in self.newdata.keys():
            dsout.variables[v][:] = self.newdata[v][:]
        dsout.close()

    def delete_variable(self, varname):
        del self.variables[varname]
        return(self)

    def insert_variable(self, var_dict, data):
        '''
        <var_dict> is a dictionary as in self.variables.
        <data> is a dictionary {<varname>: numpy.array(...)}
        '''
        self.variables.update(var_dict)
        self.newdata.update(data)
        return(self)

    def insert_dimensions(self, dimensions):
        '''<dimensions> is an OrderedDictionary {<dimname>: <dimsize>, ...}.'''
        self.dims.update(dimensions)
        return(self)

    def modify_variable_meta(self, varname, newdtype=None,
                             newdims=None, **newattributes):
        '''
        varname (str): name of variable to modify.
        newdtype (numpy.dtype): new datatype, if applicable.
        newdims (OrderedDict): new dimensions as {dimname: size, ...},
                               if applicable.
        In case newdims are given, these dimensions are created if not
        already present, and the data will be set to a properly sized
        array filled with _FillValue.
        '''
        self.variables[varname]['attributes'].update(newattributes)
        if newdtype:
            assert(type(newdtype) == np.dtype)
            self.variables[varname]['dtype'] = newdtype
        if newdims:
            assert(type(newdims) == OrderedDict)
            missdims = set(newdims) - set(self.dims)
            self.dims.update([(d, newdims[d]) for d in missdims])
            newdimnames = tuple(newdims.keys())
            self.variables[varname]['dimensions'] = newdimnames
            newdimsizes = tuple(newdims.values())
            self.newdata.update(
                self._mk_empty_data(varname, newdimsizes,
                                    self.variables[varname]['dtype']))
        return(self)

    def modify_variable_data(self, newdata):
        '''
        (dict) newdata: new data as {<varname>: numpy.array, ...}
        Attaches <newdata> to <varname>.
        '''
        v_undef = list(set(newdata.keys()) - set(self.variables.keys()))
        v_def = list(set(newdata.keys()) & set(self.variables.keys()))       
        if v_undef:
            print("WARNING: data attached to non-existing variables {}"
                  .format(v_undef))
        if v_def:
            # set unlimited dimensions to None
            shapes_expect = [(varname,
                              self._get_dimshape(varname),
                              newdata[varname].shape,
                              self.variables[varname]['dtype'],
                              newdata[varname].dtype) for varname in v_def]
            print(shapes_expect)
#############################################
# compare only dimension elements where shapes_expect[1] has not a None
            # mismatch = [x[0] for x in shapes_expect if (x[1] is not None and
            #                                             x[1] != x[2])]
# ##########################################
            if mismatch:
                print("WARNING: Dimension mismatch for variables: {}."
                      .format(mismatch))
            mismatch = [x[0] for x in shapes_expect if x[3] != x[4]]
            if mismatch:
                print("WARNING: Datatype mismatch for variables: {}"
                      .format(mismatch))
        self.newdata.update(newdata)
        return(self)


class Compress(NcFilter):
    def __init__(self, origin):
        super(Compress, self).__init__(origin)
        self.outResolutionShort = 2.0**16 - 2
        self.outResolutionLong = 2.0**32 - 2

    def _compress_prep(self, vname):
        '''
        Prepare lossy compression of variable <vname>.
        Check range, computed offset and scaling, and check if variable is
        well behaved (short integer ok) or highly skewed (long integer necessary).
        Return parameters for compressing.
        '''
        v = self._get_origin_values(vname)
        minVal = np.min(v[:])
        maxVal = np.max(v[:])
        meanVal = np.mean(v[:])
        if np.min([meanVal - minVal,
                   maxVal - meanVal]) < (maxVal - minVal) / 1000.:
            intType = np.dtype('uint32')
            outres = self.outResolutionLong
            fillval = np.uint32(2**32 - 1)
        else:
            intType = np.dtype('uint16')
            outres = self.outResolutionShort
            fillval = np.uint16(2**16 - 1)
        return(minVal, meanVal, maxVal, v.dtype, intType, outres, fillval)

    def _update_history_att(self):
        newhistory = (datetime.datetime.now().ctime() +
                      ': ' + ' '.join(sys.argv))
        try:
            newatt = "{}\n{}".format(newhistory, self.glob_atts['history'])
            #  separating new entries with "\n" because there is an undocumented
            #  feature in ncdump that will make it look like the attribute is an
            #  array of strings, when in fact it is not.
        except AttributeError:
            newatt = newhistory
        self.glob_atts['history'] = newatt



    


    # def compress(self):
    #     

    #         def _select_vars():
    #             '''Select variables that are going to be packed'''
    #             v_sel = [x for x in self.ds.variables.iteritems()
    #                      if (x[0] not in P.exclude) and
    #                      (x[1].ndim >= 2) and
    #                      (x[1].dtype in ['float64', 'float32', 'uint32', 'uint16'])]
    #             self.selected_vars = v_sel




            
        
        
                                
    #                  for x in self.dsin.variables.values()]
    #     for v in self.dsin.variables.itervalues():
    #         v_new = self.dsout.createVariable(v.name, v.dtype, v.dimensions)
    #             atts_new = dict([(x, v.getncattr(x)) for x in v.ncattrs()])
    #             v_new.setncatts(atts_new)
    #         else:
    #             v_new = compresshook(v)
    #         v_new[:] = v[:]
        
        
        

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
       # self.fin, self.fout, self.overwrite = self.parse_cmd()
       # self.dsout = Dataset(self.fout, 'w')





    def cp_all(self, compressvars=None, compresshook=None):
        '''
        Copy content of netCDF-structure from self.dsin to self.dsout. Replace
        variables in <compressvars> with the output of compresshook(<variable>).
        '''
        # Global attributes
        glob_atts = dict([(x, self.dsin.getncattr(x))
                          for x in self.dsin.ncattrs()])
        self.dsout.setncatts(glob_atts)
        # dimensions
        dim_sizes = [None if x.isunlimited() else len(x)
                     for x in self.dsin.dimensions.values()]
        dimensions = zip(self.dsin.dimensions.keys(), dim_sizes)
        for d in dimensions:
            self.dsout.createDimension(d[0], size=d[1])
        # variables
        for v in self.dsin.variables.itervalues():
            print("processing variable: {}".format(v.name)),
            if compressvars is None or v.name not in compressvars:
                print("copy")
                v_new = self.dsout.createVariable(v.name, v.dtype, v.dimensions)
                atts_new = dict([(x, v.getncattr(x)) for x in v.ncattrs()])
                v_new.setncatts(atts_new)
            else:
                v_new = compresshook(v)
            v_new[:] = v[:]

    def check_values(self, v):
        '''
        Checks whether values of <v> are identical for self.dsin
        and self.dsout.
        '''
        assert(np.all(self.dsin.variables[v][:] == self.dsout.variables[v][:]))
        return("Value check for {} passed.".format(v))

    # def compress(self, v):

    # outResolutionShort = 2.0**16 - 2
    # outResolutionLong = 2.0**32 - 2  # for unknown reason 2**32 produces wrong results
    #                                  # try it anyway - hvw
    # complevel = 9

    # # coordinate variables to prevent from compression even if they are 2D
    # # TODO automatically find coordinate variables!
    # exclude = ('lon', 'lat', 'slon', 'slat', 'slonu', 'slatu', 'slonv',
    #            'slatv', 'time', 'time_bnds', 'rlon', 'rlat', 'level_bnds',
    #            'level', 'levels')

    #     # check range, computed offset and scaling, and check if variable is
    #     # well behaved (short integer ok) or highly skewed (long integer necessary)
    #     minVal = np.min(v[:])
    #     maxVal = np.max(v[:])
    #     meanVal = np.mean(v[:])
    #     if np.min(meanVal - minVal,
    #               maxVal - meanVal) < (maxVal - minVal) / 1000:
    #         intType = np.dtype('uint32')
    #         outres = self.outResolutionLong
    #         fillval = np.uint32(2**32 - 1)
    #     else:
    #         intType = np.dtype('uint16')
    #         outres = self.outResolutionShort
    #         fillval = np.uint16(2**16 - 1)
    #     print("Packing variable {} [min:{}, mean:{}, max:{}] <{}> into <{}>"
    #           .format(v.name, minVal, meanVal, maxVal, v.dtype, intType))

    #     # choose chunksize: The horizontal domain (last 2 dimensions)
    #     # is one chunk. That the last 2 dimensions span the horizontal
    #     # domain is a COARDS convention, which we assume here nonetheless.
    #     chunksizes = tuple([1]*(len(v.dimensions) - 2) +
    #                        [len(self.dsin.dimensions[x])
    #                         for x in v.dimensions[-2:]])
    #     v_new = self.dsout.createVariable(v.name, intType, v.dimensions,
    #                                       zlib=True, complevel=self.complevel,
    #                                       chunksizes=chunksizes,
    #                                       fill_value=fillval)
    #     scale_factor = (maxVal - minVal) / outres or 1
    #     v_new.setncattr('scale_factor', scale_factor)
    #     v_new.setncattr('add_offset', minVal)
    #     v_new.setncattr('_FillValue', fillval)
    #     v_new.set_auto_maskandscale(True)
    #     # copy untouched attributes
    #     att_cp = dict([(x, v.getncattr(x)) for x in v.ncattrs()
    #                    if x not in v_new.ncattrs()])
    #     v_new.setncatts(att_cp)
    #     return(v_new)

    def get_coordvars(self, dimension=None, type=None):
        '''
        Returns all variable names that represent coordinates. Restrict
        to time, easting, northing by specifying
          dimension='T',
          dimension='Z',
          dimension='X',
          dimension='Y',
        respectively.
        Specify
          <type>='dim'
        to get only dimension variables, or
          <type>='aux'
        to get only auxiliary coordinate variables.
        '''
        def isT(v):
            print(v.name)
            try:
                if v.getncattr('axis') == 'T':
                    return(True)
            except:
                pass
            try:
                if (v.getncattr('units').split(' ')[0]
                    in ['common_year', 'common_years', 'year', 'years', 'yr', 'a', 'month', 'months', 'week',
                        'weeks', 'day', 'days', 'd', 'hour', 'hours', 'hr',
                        'h', 'minute', 'minutes', 'min', 'second', 'seconds',
                        's', 'sec']):
                    return(True)
            except:
                pass
            return(False)

        


if __name__ == "__main__":
    P = NcFilter('test_unpacked.nc')
    P.write('test_packed.nc')
    #P.dsout.close()
    # print (P.fin, P.fout, P.overwrite)
    # P.cp_all(compressvars='pr', compresshook=P.compress)
    # # print(P.check_values('pr'))
    # P.dsout.close()
    # P.dsin.close()
    

# run compress_netcdf.py test_unpacked.nc -o test_packed.nc
