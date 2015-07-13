from compress_netcdf import *
from nose.tools import *
import os
import gzip
from urllib2 import urlopen

TESTIN = 'DMI-HIRHAM5_A1B_ARPEGE_MM_25km_pr.nc'
TESTOUT = 'testout.nc'


class NcFilter_Test():

    def setUp(self):
        try:
            os.remove(TESTOUT)
        except OSError:
            pass
        assert not os.path.exists(TESTOUT)
        if not os.path.exists(TESTIN):
            blk = 2**20
            print('Downloading testfile ...'),
            f_rem = urlopen('http://ensemblesrt3.dmi.dk/data/A1B/DMI/MM/'
                            + TESTIN + '.gz')
            with open(TESTIN + '.gz', 'w') as fout:
                while True:
                    block = f_rem.read(blk)
                    if not block:
                        break
                    fout.write(block)
            print(' finished!')
            with gzip.open(TESTIN + '.gz', 'r') as f_in:
                f_out = open(TESTIN, 'w')
                while True:
                    block = f_in.read(blk)
                    if not block:
                        break
                    f_out.write(block)
                f_out.close()
        self.P = NcFilter(TESTIN)

    def write_meta_test(self):
        self.P.write(TESTOUT)
        assert(self._comparemeta(TESTOUT))
        self.P.glob_atts['fake'] = 1000
        assert(not self._comparemeta(TESTOUT))

    def _comparemeta(self, file2):
        P2 = NcFilter(file2)
        return(self.P.glob_atts == P2.glob_atts and
               self.P.dims == P2.dims and
               self.P.variables == P2.variables)

    def write_data_cp_test(self):
        self.P.write(TESTOUT)
        ds1 = Dataset(TESTIN, 'r')
        ds2 = Dataset(TESTOUT, 'r')
        for v in ds1.variables:
            assert (np.all(ds1.variables[v][:] == ds2.variables[v][:]))
        ds1.close()
        ds2.close()

    def write_data_newvar_test(self):
        self.P.variables.update({'newvar': {
            'dtype': 'int',
            'dimensions': ('newdim1',
                           'newdim2',
                           'newdim3'),
            'attributes': {}}})
        self.P.dims['newdim1'] = 2
        self.P.dims['newdim2'] = 3
        self.P.dims['newdim3'] = 4
        self.P.newdata = {'newvar': np.arange(1, 25).reshape(2, 3, 4)}
        self.P.write(TESTOUT)
        self._comparemeta(TESTOUT)
        ds2 = Dataset(TESTOUT, 'r')
        assert(np.all(self.P.newdata['newvar'] == ds2.variables['newvar'][:]))
        ds2.close()

    def delete_variable_test(self):
        self.P.delete_variable('lon').write(TESTOUT)
        ds1 = Dataset(TESTIN, 'r')
        ds2 = Dataset(TESTOUT, 'r')
        assert((set(ds1.variables.keys()) - set(ds2.variables.keys()))
               == {'lon'})

    def insert_variable_test(self):
        newdims = {'newdim1': 4, 'newdim2': 5, 'newdim3': 6}
        self.P.insert_dimensions(newdims)
        var_dict = {'testinsert': {
            'dtype': float,
            'dimensions': ('newdim1', 'newdim2', 'newdim3'),
            'attributes': {'att1': 1, 'att2': 'two', 'att3': 3.01}}
        }
        data = {'testinsert': np.random.randn(4, 5, 6)}
        self.P.insert_variable(var_dict, data).write(TESTOUT)
        d1 = Dataset(TESTOUT, 'r')
        d1v = d1.variables['testinsert']
        assert(np.all(d1v[:] == data['testinsert']))
        d1.close()
        assert(self._comparemeta(TESTOUT))

    def modify_variable_meta_test(self):
        # '''
        # newattributes: old not mentioned -> keep old,
        # old = None -> delete
        # old = value -> replace
        # new = value -> insert)
        # '''
        newdims = OrderedDict([('newdim1', 4), ('newdim2', 5), ('newdim3', 6)])
        newdimensions = tuple(newdims.keys())
        newdimshape = tuple(newdims.values())
        @raises(AssertionError)
        def fail_dimensions_test():
            self.P.modify_variable_meta('pr', newdims=newdimensions,
                                        units='buckets per squarefoot',
                                        new_att='newatt').write(TESTOUT)
        fail_dimensions_test()
        self.P.modify_variable_meta('pr', newdims=newdims,
                                    units='buckets per squarefoot',
                                    new_att='newatt').write(TESTOUT)
        d1 = Dataset(TESTOUT, 'r')
        assert(d1.variables['pr'].dimensions == newdimensions)
        assert(d1.variables['pr'].getncattr('units') == 'buckets per squarefoot')
        assert(d1.variables['pr'].getncattr('new_att') == 'newatt')
        print(d1.variables['pr'][:].shape)
        assert(d1.variables['pr'][:].shape == newdimshape)
        d1.close()
    
    def modify_variable_meta_dtype_test(self):
        self.P.modify_variable_meta('pr', newdtype=np.dtype('uint16'),
                                    _FillValue=None, missing_value=None)
        self.P.write(TESTOUT)
        d1 = Dataset(TESTOUT, 'r')
        assert(d1.variables['pr'].dtype == np.dtype('uint16'))

    def modify_variable_data_test(self):
        newdata = {'rlat': np.arange(190, dtype='float32')}  # OK
        newdata.update({'xx': np.arange(5)})  # non-exist
        newdata.update({'lat': np.zeros((3, 3))})  # wrong shape
        self.P.modify_variable_data(newdata)
        output = sys.stdout.getvalue().strip()
        assert("WARNING: data attached to non-existing variables ['xx']"
               in output)
        assert("WARNING: Dimension mismatch for variables: ['lat']" in output)
        assert("WARNING: Datatype mismatch for variables: ['lat']" in output)
        assert("['rlat']" not in output)
        self.P.newdata = {}
        newdata = {'rlat': np.arange(190, dtype='float32')}
        self.P.modify_variable_data(newdata).write(TESTOUT)
        assert(np.all(Dataset(TESTOUT, 'r').variables['rlat'][:]
                      == np.arange(190, dtype='float32')))

    def insert_dimensions_test(self):
        newdims = {'newdim1': 4, 'newdim2': 5, 'newdim3': 6}
        self.P.insert_dimensions(newdims).write(TESTOUT)
        P2 = NcFilter(TESTOUT)
        assert(self.P.dims == P2.dims)

    def _get_dimshape_test(self):
        ds1 = self.P._get_dimshape('pr')
        ds2 = self.P._get_dimshape('time')
        ds3 = self.P._get_dimshape('rlon')
        print(ds1, ds2, ds3)
        assert(ds1 == (None, 190, 174) and ds2 == (None, )
               and ds3 == (174, ))

    def _get_origin_values_test(self):
        pr = self.P._get_origin_values('pr')
        pr1 = Dataset(TESTIN, 'r').variables['pr'][:]
        assert(np.all(pr == pr1))

    def update_history_att_test(self):
        self.P.glob_atts['history'] = "oldhistory attribute"
        newhistory = (datetime.datetime.now().ctime() +
                      ': ' + ' '.join(sys.argv))
        self.P.update_history_att()
        assert(self.P.glob_atts['history']
               == "{}\n{}".format(newhistory, "oldhistory attribute"))
        del self.P.glob_atts['history']
        newhistory = (datetime.datetime.now().ctime() +
                      ': ' + ' '.join(sys.argv))
        self.P.update_history_att()
        assert(self.P.glob_atts['history'] == newhistory)


import scipy.stats as scst


class Compress_Test():
    def setUp(self):
        try:
            os.remove(TESTOUT)
        except OSError:
            pass
        assert not os.path.exists(TESTOUT)
        self.C = Compress(TESTIN)

    def _compress_prep_small_test(self):
        ret = self.C._compress_prep('pr')
        v1 = self.C._get_origin_values('pr')
        des = scst.describe(v1, axis=None)
        # print(ret)
        # print(des)
        assert(ret[0:3] == (des.minmax[0], des.mean, des.minmax[1]))
        assert(ret[4:7] == (2.0**16 - 2, np.dtype('uint16'),
                            np.uint16(2**16 - 1)))
        assert(ret[3] == ((des.minmax[1] - des.minmax[0]) / 2.0**16 - 2) or 1)

    def _compress_prep_big_test(self):
        v1 = self.C._get_origin_values('pr')
        des = scst.describe(v1, axis=None)
        repmax = 1000 * des.mean - 999 * des.minmax[0] + 1
        v1.flat[np.argmax(v1)] = repmax
        newdata = {'pr': v1}
        self.C.modify_variable_data(newdata).write(TESTOUT)
        C1 = Compress(TESTOUT)
        ret = C1._compress_prep('pr')
        v1 = C1._get_origin_values('pr')
        des = scst.describe(v1, axis=None)
        assert(ret[0:3] == (des.minmax[0], des.mean, des.minmax[1]))
        assert(ret[4:7] == (2.0**32 - 2, np.dtype('uint32'),
                            np.uint32(2**32 - 1)))
        assert(ret[3] == ((des.minmax[1] - des.minmax[0]) / 2.0**32 - 2) or 1)

    def _find_compressible_variables_test(self):
        compvars, excludevars = self.C._find_compressible_variables()
        assert(compvars == [u'pr'])
        assert(excludevars == [u'rlat', u'rlon', u'rotated_pole', u'time',
                               u'lon', u'lat', u'time_bnds', 'slon', 'slat',
                               'slonu', 'slatu', 'slonv', 'slatv',
                               'level_bnds', 'level', 'levels'])

    def _calc_chunksizes_test(self):
        res = self.C. _calc_chunksizes('pr')
        assert(res == [1, 190, 174])

    def compress_test(self):
        cparams = self.C._compress_prep('pr')
        maxerr = np.abs(0.5 * cparams[3])
        print("maxerr: {}".format(maxerr))
        self.C.compress().write(TESTOUT)
        dout = Dataset(TESTOUT, 'r').variables['pr'][0:100, :]
        din = Dataset(TESTIN, 'r').variables['pr'][0:100, :]
        maxdiff = np.max(abs(dout - din))
        print("maxdiff: {}".format(maxdiff))
        assert(maxdiff <= maxerr)
        
        # self.C.compress()
        #print("\n returned: {}".format(compvars))
        # raise Exception
        
        

    # def _compress_prep_test(self):
    #     d_skewed =

        
    # def copy_variable_meta_test(self):#, varname, newname):
    #     raise Exception

    # def delete_dimensions_test(self):#, dimnames):
    #     raise Exception

    # def compress_test(self):
    #     raise Exception















