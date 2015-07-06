from compress_netcdf import *
import os

TESTIN = 'test_unpacked_w.nc'
TESTOUT = 'test_packed.nc'


class NcFilter_Test():

    def setup(self):
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

        # def write_data_replace_vars_test(self):
            
        # def write_data_remove_vars_test(self):










