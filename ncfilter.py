#!/usr/bin/env python
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
            # 'flags' is an OrderedDict with netCDF4.Variable methods as keys
            # and a list of arguments.
            # E.g. {Datset.Variable.set_auto_mask: [True]}
            self.variables = OrderedDict([(x.name, {
                'dtype': x.dtype,
                'dimensions': x.dimensions,  # tuple
                'attributes': self._get_var_attrs(x),
                'flags': OrderedDict(),
                'createargs': OrderedDict()})
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

    def write(self, outfile, histstring="_undef_"):
        '''
        Creates <outfile> with meta-data as in
          self.glob_atts
          self.dims
          self.variables
        Data in self.newdata ( {<varname>: np.array(...), ...} ) replaces
        the respective data of <varname> in the original file.
        New <varname>s in <self.newdata> have to be present in <self.variables>.

        <histstring> will be prepended to the global attribute "history".
        If no such attribute exists, it will be created.
        Set histstring=None to leave the attribute untouched. Not recommended!
        '''
        renam = False
        if self.origin == outfile:
            origout = outfile
            outfile += '.tmp'
            renam = True
        dsout = Dataset(outfile, "w")
        dsout.set_auto_mask(True)
        
        # add to history attribute
        self.update_history_att(histstring)

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
                                        dimensions=v['dimensions'],
                                        **v.get('createargs', {}))
            vout.setncatts(v['attributes'])
            for f in v.get('flags', {}).iteritems():
                getattr(vout, f[0])(*f[1])

        # variables to be identically copied (data):
        vcp = set(self.variables.keys()) - set(self.newdata.keys())
        with Dataset(self.origin, "r") as dsin:
            for v in vcp:
                dsout.variables[v][:] = dsin.variables[v][:]

        # variables with new data
        for v in self.newdata.keys():
            dsout.variables[v][:] = self.newdata[v][:]
        dsout.close()
        # replace original file with temporary output
        if renam:
            os.rename(outfile, origout)

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

        The remaining named parameters **newattributes update (i.e. append
        and / overwrite the attributes of <varname>. A named parameter wit
        value = None will result in the deletion of the attribute.
        '''
        self.variables[varname]['attributes'].update(newattributes)
        for atnam, atval in self.variables[varname]['attributes'].items():
            if atval is None:
                del self.variables[varname]['attributes'][atnam]
        if newdtype is not None:
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
            mismatch = {}
            for m in shapes_expect:
                if (len(m[1]) != len(m[2]) or
                    not np.all([x == y or None in [x, y]
                                for x, y in zip(m[1], m[2])])):
                    mismatch[m[0]] = (
                        "WARNING: dimensions don't match: {} vs. {}"
                        .format(m[1], m[2]))
            if mismatch:
                print(mismatch)
                print("Shapes expect: {}".format(shapes_expect))
            mismatch = [x[0] for x in shapes_expect if x[3] != x[4]]
            if mismatch:
                print("WARNING: Datatype mismatch for variables: {}"
                      .format(mismatch))
        self.newdata.update(newdata)
        return(self)

    def update_history_att(self, newhist="_undef_"):
        '''Precedes current global attribute "history" with date + command'''

        if __name__ == "__main__":
            newhistory = (datetime.datetime.now().ctime() +
                          ': ' + ' '.join(sys.argv))
        elif newhist == "_undef_":
            print("Warning: No new history attribute given. Using 'unspecified action'")
            newhistory = (datetime.datetime.now().ctime() +
                          ': ' + "unspecified action")
        elif newhist is None:
            print("Warning: History attribute left unchanged!")
            return(self)
        else:
            newhistory = newhist
        try:
            newatt = "{}\n{}".format(newhistory, self.glob_atts['history'])
            #  separating new entries with "\n" because there is an undocumented
            #  feature in ncdump that will make it look like the attribute is an
            #  array of strings, when in fact it is not.
        except KeyError:
            newatt = newhistory
        self.glob_atts['history'] = newatt
        return(self)

    def checkarg(self):
        print("\nsys.argv: {}".format(sys.argv))
        print("__name__: {}".format(__name__))


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
        # scale factor = 1 if maxVal == minVal
        scale_factor = (maxVal - minVal) / outres or 1
        return(minVal, meanVal, maxVal, scale_factor, outres, intType, fillval)

    def _find_compressible_variables(self):
        ''' Returns variable names that are not thought
        to be coordinate variables'''

        # It is quite difficult to properly identify the coordinate variables
        # assuming CF-Conventions (1.6) only. Therefore assume all 1-D variables
        # need not be compressed.
        # exclude proper coordinate variables (1-dimensional)
        exclude_coord = [varname for varname in self.variables if
                         len(self.variables[varname]['dimensions']) <= 1]
        # exclude auxiliary coordinates and cell-bounds
        exclude_aux_coords = []
        for atts in [v['attributes'] for v in self.variables.values()]:
            auxcoords = atts.get('coordinates') or ''
            auxcoords += ' ' + (atts.get('bounds') or '')
            exclude_aux_coords.extend(auxcoords.split())
        # for good measure exclude variable names from Dominik's list
        exclude_dom = ['lon', 'lat', 'slon', 'slat', 'slonu', 'slatu', 'slonv',
                       'slatv', 'time', 'time_bnds', 'rlon', 'rlat',
                       'level_bnds', 'level', 'levels']
        # also exclude variables of wrong datatype
        exclude_dtyp = []
        comp_dtyp = [np.dtype(x) for x in ['float64', 'float32',
                                           'uint32', 'uint16']]
        for vn, v in self.variables.iteritems():
            if v['dtype'] not in comp_dtyp:
                exclude_dtyp.append(vn)
        exclude_all = exclude_coord + exclude_aux_coords + \
                      exclude_dom + exclude_dtyp
        exclude_all = list(OrderedDict.fromkeys(exclude_all))  # make unique
        compressible = [v for v in self.variables if v not in exclude_all]
        return((compressible, exclude_all))

    def _calc_chunksizes(self, varname):
        # choose chunksize: The horizontal domain (last 2 dimensions)
        # is one chunk. That the last 2 dimensions span the horizontal
        # domain is a COARDS convention, which we assume here nonetheless.
        chu0 = [1] * (len(self.variables[varname]['dimensions']) - 2)
        chu1 = [self.dims[x] for x in self.variables[varname]['dimensions'][-2:]]
        chunksizes = chu0 + chu1
        return(chunksizes)

    def compress(self, complevel=9):
        for varname in self._find_compressible_variables()[0]:
            minVal, meanVal, maxVal,\
                scale_factor, outres, intType, fillval = self._compress_prep(varname)
            # cast fillval to new integer type
            fillval = np.array([fillval], dtype=intType)[0]
            chunksizes = self._calc_chunksizes(varname)
            # set new dType, set(reset) appropriate attributes
            self.modify_variable_meta(varname, newdtype=intType,
                                      scale_factor=scale_factor,
                                      add_offset=minVal,
                                      _FillValue=fillval)
            if 'missing_value' in self.variables[varname]['attributes']:
                self.modify_variable_meta(varname, missing_value=fillval)
            # setting 'set_suto_scale' to False, (and the packing is done
            # explicitely, because the automatic
            # packing truncates to int instead of rounding.
            self.variables[varname]['flags'] = OrderedDict([
                ('set_auto_mask', [True]),
                ('set_auto_scale', [False])])
            # set parameters for netCDF4.Dataset.createVariable 
            # (zlib-compression and fill-vale)
            self.variables[varname]['createargs'] = OrderedDict([
                ('zlib', True), ('complevel', complevel),
                ('chunksizes', chunksizes), ('fill_value', fillval)])
            # pack variable
            newdata = np.round(
                (self._get_origin_values(varname) - minVal)
                / scale_factor).astype(intType)
            self.modify_variable_data({varname: newdata})
        return(self)


class Commands(object):

    @staticmethod
    def delvar(argparse):
        argdict = argparse[0]
        parser = argparse[1]
        try:
            delvar = argdict['arguments'][0]
        except:
            parser.error("delvar requires one ARG that names the variable to delete.")
        NcFilter(argdict['fin']).delete_variable(delvar).write(argdict['fout'])

    @staticmethod
    def compress(argparse):
        argdict = argparse[0]
        parser = argparse[1]
        try:
            cl = int(argdict['arguments'][0])
        except:
            cl = 9
        Compress(argdict['fin']).compress(complevel=cl).write(argdict['fout'])


def main():
    def _get_commands():
        return([m for m in dir(Commands) if not m.startswith('__')])

    parser = argparse.ArgumentParser(description='Performs operations' +
                                     ' on a netCDF file.',
                                     epilog='OUTFILE = INFILE will replace the ' +
                                     'INFILE with the OUTFILE.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('command', help="possible commands: {}"
                        .format(_get_commands()), metavar='COMMAND')
    parser.add_argument('arguments', help='arguments for commands' +
                        '\n    "compress": ARG = compression_level (int), default=9' +
                        '\n    "delvar": ARG = variable_to_delete (str)',
                        metavar="ARG", nargs='*')
    parser.add_argument('fin', help='input file', metavar='INFILE')
    parser.add_argument('fout', help='output file',
                        metavar='OUTFILE')
    args = vars(parser.parse_args())

    # check input file
    if not os.access(args['fin'], os.R_OK):
        parser.error("Can't open {} for reading".format(args['fin']))

    # check output file
    outpath = os.path.dirname(args['fin']) or '.'
    if not os.access(outpath, os.W_OK):
        parser.error("can't write output file {}".format(args['fout']))

    # check command
    if not hasattr(Commands, args['command']):
        parser.error("Command {} not implemented".format(args['command']))
    else:
        getattr(Commands, args['command'])([args, parser])

if __name__ == "__main__":
    main()
