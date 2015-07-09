from compress_netcdf import *
from nose.tools import *
import os

TESTIN = 'test_unpacked_w.nc'
TESTOUT = 'test_packed.nc'


class NcFilter_Test():

    def setUp(self):
        try:
            os.remove(TESTOUT)
        except OSError:
            pass
        assert not os.path.exists(TESTOUT)
        self.P = NcFilter(TESTIN)

    def no_test(self):
        pass

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
        self.P.variables.append({'name': 'newvar',
                                 'dtype': 'int',
                                 'dimensions': ('newdim1',
                                                'newdim2',
                                                'newdim3'),
                                 'attributes': {}})
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
        var_dict = {'name': 'testinsert', 'dtype': float,
                    'dimensions': ('newdim1', 'newdim2', 'newdim3'),
                    'attributes': {'att1': 1, 'att2': 'two', 'att3': 3.01}}
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
    # def copy_variable_meta_test(self):#, varname, newname):
    #     raise Exception

    # def delete_dimensions_test(self):#, dimnames):
    #     raise Exception

    # def compress_test(self):
    #     raise Exception






