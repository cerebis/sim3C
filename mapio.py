import numpy as np
import pandas
import h5py
import io_utils


def write_map(data, fname, fmt, names=None):
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
    :param names: names used for both rows an columns
    """

    # try to convert sparse matrices to a dense form before writing
    if not isinstance(data, np.ndarray):
        try:
            data = data.todense()
        except AttributeError:
            raise RuntimeError('Supplied data object does not appear to be a sparse matrix, nor a numpy array')

    if fmt == 'csv':
        # append a suffix is not found.
        if not fname.endswith('csv'):
            fname = '{}.csv'.format(fname)
        if names is None:
            # plain CSV, no row or column names
            pandas.DataFrame(data).to_csv(fname, header=False, index=False, sep=',')
        else:
            # add row and column names to CSV if supplied
            pandas.DataFrame(data, index=names, columns=names).to_csv(fname, sep=',')

    elif fmt == 'h5':
        if not fname.endswith('h5'):
            fname = '{}.h5'.format(fname)
        with h5py.File(fname, 'w') as hndl:
            hndl.create_dataset('map', data=data, compression='gzip')
            # add names as as separate entry
            if names is not None:
                hndl.create_dataset('names', data=names, compression='gzip')

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
    data = None
    cols = None

    if fmt == 'csv':
        with io_utils.open_input(fname) as csv_input:
            try:
                if names:
                    df = pandas.read_csv(csv_input, delimiter=delim, index_col=0, header=0)
                    data = df.as_matrix().astype(np.float)
                    if (df.index != df.columns).all():
                        raise IOError('row and column names do not agree')
                    cols = np.asarray(df.columns.tolist())
                else:
                    df = pandas.read_csv(csv_input, delimiter=delim, dtype=np.float, header=None)
                    data = df.as_matrix()
            except ValueError as er:
                print 'Invalid data. Expected numeric values. Perhaps wrong delimiter used'
                raise er

    elif fmt == 'h5':
        with h5py.File(fname, 'r') as hndl:
            data = hndl['map'][:]
            if names:
                cols = np.asarray(hndl['names'][:])
    else:
        raise RuntimeError('Unimplemented format [{}]'.format(fmt))

    return data, cols

