from compress_netcdf import *
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
        newdata = {'newvar': np.arange(1, 25).reshape(2, 3, 4)}
        self.P.write(TESTOUT, newdata=newdata)
        self._comparemeta(TESTOUT)
        ds2 = Dataset(TESTOUT, 'r')
        assert(np.all(newdata['newvar'] == ds2.variables['newvar'][:]))
        ds2.close()

    def delete_variable_test(self): #, varname):
        self.P.delete_variable('lon').write(TESTOUT)
        ds1 = Dataset(TESTIN, 'r')
        ds2 = Dataset(TESTOUT, 'r')
        assert((set(ds1.variables.keys()) - set(ds2.variables.keys()))
               == {'lon'})

    def insert_variable_test(self): #, var_dict, data):
        raise Exception

    def modify_variable_meta_test(self):  #, varname, newname=None, newdtype=None, newdimensions=None, **newattributes):
        # '''
        # newattributes: old not mentioned -> keep old,
        # old = None -> delete
        # old = value -> replace
        # new = value -> insert)
        # '''
        raise Exception

    def copy_variable_meta_test(self):#, varname, newname):
        raise Exception

    def replace_variable_data_test(self):#, varname, new_data_func=None,
                                   # *newdataargs):
        raise Exception

    def delete_dimensions_test(self):#, dimnames):
        raise Exception

    def insert_dimensions_test(self):#, dimensions):
        raise Exception

    def compress_test(self):
        raise Exception






