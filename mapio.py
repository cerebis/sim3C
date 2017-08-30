import numpy as np
import pandas
import h5py


def write_map(data, fname, fmt, field_as=np.float):
    """
    Write a contact map (2D numpy array) to file. File formats
    are either CSV or a simple HFD5 with a single dataset named
    'map'. 
    
    A filename suffix is appended automatically if omitted:
    'csv' or 'h5'.
    
    In each case, currently all matrices are converted to a dense 
    representation before saving. This may require substantial memory.
    
    :param data: 2d numpy array to write
    :param fname: output file name
    :param fmt: file format 'csv' or 'h5'
    :param field_as: number representation of each element
    """
    field_fmt = '%f'
    if field_as == np.int:
        data = data.astype(np.int)
        field_fmt = '%d'

    # try to convert sparse matrices to a dense form before writing
    if not isinstance(data, np.ndarray):
        try:
            data = data.todense()
        except AttributeError:
            raise RuntimeError('Supplied data object does not appear to be a sparse matrix, nor a numpy array')

    if fmt == 'csv':
        if not fname.endswith('csv'):
            fname = '{}.csv'.format(fname)
        np.savetxt(fname, data, fmt=field_fmt, delimiter=',')

    elif fmt == 'h5':
        if not fname.endswith('h5'):
            fname = '{}.h5'.format(fname)
        with h5py.File(fname, 'w') as hndl:
            hndl.create_dataset('map', data=data, compression='gzip')

    else:
        raise RuntimeError('Unimplemented format [{}]'.format(fmt))


def read_map(fname, fmt, delim=',', names=False):
    """
    Read a contact map from file. The map is returned as a 2D
    numpy array. File formats are 'CSV' or a simple HDF5 with
    a single dataset named 'map'. In the case of a named CSV,
    a tuple will be returned (matrix, rownames, colnames).
    :param fname: input file name
    :param fmt: file format 'csv' or 'h5'
    :param delim: optional delimiter for csv
    :param names: for CSV, the first row and column specify names.
    :return: 
    """
    if fmt == 'csv':
        if names:
            df = pandas.read_csv(fname, delimiter=delim, index_col=0, header=0)
            return df.as_matrix(), df.index.tolist(), df.columns.tolist()
        else:
            data = np.loadtxt(fname, delimiter=delim)
    elif fmt == 'h5':
        if names:
            print 'Warning: names parameter ignored for h5 format'
        with h5py.File(fname, 'r') as hndl:
            data = hndl['map'][:]
    else:
        raise RuntimeError('Unimplemented format [{}]'.format(fmt))
    return data

