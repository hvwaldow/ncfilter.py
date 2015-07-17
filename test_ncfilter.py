#
# test with py.test
# http://pytest.org
#

from ncfilter import *
import os
import gzip
from urllib2 import urlopen
import pytest

TESTIN = 'DMI-HIRHAM5_A1B_ARPEGE_MM_25km_pr.nc'
TESTOUT = 'testout.nc'


@pytest.fixture(scope='module')
def get_nctestfile():
    if not os.path.exists(TESTIN):
        blk = 2**20
        print('Downloading testfile {} ... '.format(TESTIN)),
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


@pytest.fixture()
def del_ncout_testfile():
    try:
        os.remove(TESTOUT)
    except OSError:
        pass
    assert not os.path.exists(TESTOUT)


@pytest.fixture()
def ncfilter(get_nctestfile, del_ncout_testfile):
    P = NcFilter(TESTIN)
    return(P)


@pytest.fixture()
def compress(get_nctestfile, del_ncout_testfile):
    C = Compress(TESTIN)
    return(C)


def test_write_meta(ncfilter):
    ncfilter.write(TESTOUT)
    assert(_comparemeta(ncfilter, TESTOUT))
    ncfilter.glob_atts['fake'] = 1000
    assert(not _comparemeta(ncfilter, TESTOUT))


def _comparemeta(ncfilter, file2):
    P2 = NcFilter(file2)
    # add non-required atributes to P2.variables
    # which get written in any case:
    for v in ncfilter.variables.values():
        try:
            a = v['createargs']
        except KeyError:
            v['createargs'] = OrderedDict()
        try:
            a = v['flags']
        except KeyError:
                v['flags'] = OrderedDict()
    return(ncfilter.glob_atts == P2.glob_atts and
           ncfilter.dims == P2.dims and
           ncfilter.variables == P2.variables)


def test_write_data_cp(ncfilter):
    ncfilter.write(TESTOUT)
    with Dataset(TESTIN, 'r') as ds1, Dataset(TESTOUT, 'r') as ds2:
        for v in ds1.variables:
            assert (np.all(ds1.variables[v][:] == ds2.variables[v][:]))


def test_write_data_newvar(ncfilter):
    ncfilter.variables.update({'newvar': {
        'dtype': 'int',
        'dimensions': ('newdim1',
                       'newdim2',
                       'newdim3'),
        'attributes': {}}})
    ncfilter.dims['newdim1'] = 2
    ncfilter.dims['newdim2'] = 3
    ncfilter.dims['newdim3'] = 4
    ncfilter.newdata = {'newvar': np.arange(1, 25).reshape(2, 3, 4)}
    ncfilter.write(TESTOUT)
    _comparemeta(ncfilter, TESTOUT)
    with Dataset(TESTOUT, 'r') as ds2:
        assert(np.all(ncfilter.newdata['newvar']
                      == ds2.variables['newvar'][:]))


def test_delete_variable(ncfilter):
    ncfilter.delete_variable('lon').write(TESTOUT)
    with Dataset(TESTIN, 'r') as ds1, Dataset(TESTOUT, 'r') as ds2:
        assert((set(ds1.variables.keys()) - set(ds2.variables.keys()))
               == {'lon'})


def test_insert_variable(ncfilter):
    newdims = {'newdim1': 4, 'newdim2': 5, 'newdim3': 6}
    ncfilter.insert_dimensions(newdims)
    var_dict = {'testinsert': {
        'dtype': float,
        'dimensions': ('newdim1', 'newdim2', 'newdim3'),
        'attributes': {'att1': 1, 'att2': 'two', 'att3': 3.01}}
    }
    data = {'testinsert': np.random.randn(4, 5, 6)}
    ncfilter.insert_variable(var_dict, data).write(TESTOUT)
    with Dataset(TESTOUT, 'r') as d1:
        d1v = d1.variables['testinsert']
        assert(np.all(d1v[:] == data['testinsert']))
    assert(_comparemeta(ncfilter, TESTOUT))


def test_modify_variable_meta(ncfilter):
    # '''
    # newattributes: old not mentioned -> keep old,
    # old = None -> delete
    # old = value -> replace
    # new = value -> insert)
    # '''
    newdims = OrderedDict([('newdim1', 4), ('newdim2', 5), ('newdim3', 6)])
    newdimensions = tuple(newdims.keys())
    newdimshape = tuple(newdims.values())
    with pytest.raises(AssertionError):
        ncfilter.modify_variable_meta('pr', newdims=newdimensions,
                                      units='buckets per squareft',
                                      new_att='newatt').write(TESTOUT)
    ncfilter.modify_variable_meta('pr', newdims=newdims,
                                  units='buckets per squareft',
                                  new_att='newatt').write(TESTOUT)
    with Dataset(TESTOUT, 'r') as d1:
        assert(d1.variables['pr'].dimensions == newdimensions)
        assert(d1.variables['pr'].getncattr('units') == 'buckets per squareft')
        assert(d1.variables['pr'].getncattr('new_att') == 'newatt')
        assert(d1.variables['pr'][:].shape == newdimshape)


def test_modify_variable_meta_dtype(ncfilter):
    ncfilter.modify_variable_meta('pr', newdtype=np.dtype('uint16'),
                                  _FillValue=None, missing_value=None)
    ncfilter.write(TESTOUT)
    with Dataset(TESTOUT, 'r') as d1:
        assert(d1.variables['pr'].dtype == np.dtype('uint16'))


def test_modify_variable_data(ncfilter, capsys):
    newdata = {'rlat': np.arange(190, dtype='float32')}  # OK
    newdata.update({'xx': np.arange(5)})  # non-exist
    newdata.update({'lat': np.zeros((3, 3))})  # wrong shape
    ncfilter.modify_variable_data(newdata)
    print(capsys.readouterr()[0].strip())
    output = capsys.readouterr()[0].strip()
    assert("WARNING: data attached to non-existing variables ['xx']"
           in output)
    assert("'lat': \"WARNING: dimensions " +
           "don't match: (190, 174) vs. (3, 3)\"" in output)
    assert("WARNING: Datatype mismatch for variables: ['lat']" in output)
    assert("['rlat']" not in output)
    ncfilter.newdata = {}
    newdata = {'rlat': np.arange(190, dtype='float32')}
    ncfilter.modify_variable_data(newdata).write(TESTOUT)
    assert(np.all(Dataset(TESTOUT, 'r').variables['rlat'][:]
                  == np.arange(190, dtype='float32')))


def test_insert_dimensions(ncfilter):
    newdims = {'newdim1': 4, 'newdim2': 5, 'newdim3': 6}
    ncfilter.insert_dimensions(newdims).write(TESTOUT)
    P2 = NcFilter(TESTOUT)
    assert(ncfilter.dims == P2.dims)


def test__get_dimshape(ncfilter):
    ds1 = ncfilter._get_dimshape('pr')
    ds2 = ncfilter._get_dimshape('time')
    ds3 = ncfilter._get_dimshape('rlon')
    assert(ds1 == (None, 190, 174) and ds2 == (None, )
           and ds3 == (174, ))


def test__get_origin_values(ncfilter):
    pr = ncfilter._get_origin_values('pr')
    pr1 = Dataset(TESTIN, 'r').variables['pr'][:]
    assert(np.all(pr == pr1))


def test_update_history_att(ncfilter, capsys):
    ncfilter.glob_atts['history'] = "oldhistory attribute"
    # no history string given
    ncfilter.update_history_att()
    assert(capsys.readouterr()[0]
           == "Warning: No new history attribute given. " +
           "Using 'unspecified action'\n")
    assert("unspecified action" in ncfilter.glob_atts["history"])
    ncfilter.glob_atts['history'] = "oldhistory attribute"
    # history string == None
    ncfilter.update_history_att(newhist=None)
    assert(capsys.readouterr()[0]
           == "Warning: History attribute left unchanged!\n")
    assert(ncfilter.glob_atts['history'] == "oldhistory attribute")


def test__compress_prep_small(compress):
    ret = compress._compress_prep('pr')
    v1 = compress._get_origin_values('pr')
    v1max, v1min, v1mean = (v1.max(), v1.min(), v1.mean())
    assert(ret[0:3] == (v1min, v1mean, v1max))
    assert(ret[4:7] == (2.0**16 - 2, np.dtype('uint16'),
                        np.uint16(2**16 - 1)))
    assert(ret[3] == ((v1max - v1min) / 2.0**16 - 2) or 1)


def test__compress_prep_big(compress):
    v1 = compress._get_origin_values('pr')
    v1max, v1min, v1mean = (v1.max(), v1.min(), v1.mean())
    repmax = 1000 * v1mean - 999 * v1min + 1
    v1.flat[np.argmax(v1)] = repmax
    newdata = {'pr': v1}
    compress.modify_variable_data(newdata).write(TESTOUT)
    C1 = Compress(TESTOUT)
    ret = C1._compress_prep('pr')
    v1 = C1._get_origin_values('pr')
    v1max, v1min, v1mean = (v1.max(), v1.min(), v1.mean())
    assert(ret[0:3] == (v1min, v1mean, v1max))
    assert(ret[4:7] == (2.0**32 - 2, np.dtype('uint32'),
                        np.uint32(2**32 - 1)))
    assert(ret[3] == ((v1max - v1min) / 2.0**32 - 2) or 1)


def test__find_compressible_variables(compress):
    compvars, excludevars = compress._find_compressible_variables()
    assert(compvars == [u'pr'])
    assert(excludevars == [u'rlat', u'rlon', u'rotated_pole', u'time',
                           u'lon', u'lat', u'time_bnds', 'slon', 'slat',
                           'slonu', 'slatu', 'slonv', 'slatv',
                           'level_bnds', 'level', 'levels'])


def test__calc_chunksizes(compress):
    res = compress. _calc_chunksizes('pr')
    assert(res == [1, 190, 174])


def test_compress(compress):
    cparams = compress._compress_prep('pr')
    compress.compress().write(TESTOUT)
    dout = Dataset(TESTOUT, 'r').variables['pr'][:]
    din = Dataset(TESTIN, 'r').variables['pr'][:]
    maxerrnorm = np.max(abs((dout - din) / cparams[3]))
    assert(maxerrnorm <= 0.51)












