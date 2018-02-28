import numpy as np
import scipy.sparse
import pandas
import h5py
import h5sparse
import io_utils


def save_h5map(fname, cmap, names=None):
    """
    Save a mixed H5 file with sparse matrix contact map and an optional list of names
    :param fname: output file name
    :param cmap: contact map in CSR, CSC or BSR format
    :param names: optional names in corresponding order to matrix elements
    """
    with h5py.File(fname, 'w') as saved:

        if cmap.format not in ('csc', 'csr', 'bsr'):
            raise NotImplementedError('Saving of sparse matrix of format {} not implemented .'.format(cmap.format))

        g = saved.create_group('map')
        g.create_dataset('indices', data=cmap.indices, compression='gzip', shuffle=True)
        g.create_dataset('indptr', data=cmap.indptr, compression='gzip', shuffle=True)
        g.create_dataset('format', data=cmap.format.encode('ascii'))
        g.create_dataset('shape', data=cmap.shape)
        g.create_dataset('data', data=cmap.data, compression='gzip', shuffle=True)

        if names is not None:
            assert len(names) == cmap.shape[0]
            saved.create_dataset('names', data=names)


def load_h5map(fname, has_names=True):
    """
    Load a mixed H5 file containing a sparse matrix and optional list of names.
    :param fname: input file name
    :param has_names: True - load the associated per-element names
    :return: cmap sparse matrix, list of names or None
    """
    with h5py.File(fname, 'r') as loaded:
        try:
            matrix_format = loaded['map/format'][()]
        except KeyError:
            raise ValueError('The file {} does not appear to contain a sparse matrix.'.format(fname))

        try:
            cls = getattr(scipy.sparse, '{}_matrix'.format(matrix_format))
        except AttributeError:
            raise ValueError('Unknown matrix format "{}"'.format(matrix_format))

        if matrix_format not in ('csc', 'csr', 'bsr'):
            raise NotImplementedError('Loading of sparse matrix of format {} not implemented.'.format(matrix_format))

        cmap = cls((loaded['map/data'][:],
                    loaded['map/indices'][:],
                    loaded['map/indptr'][:]),
                   shape=loaded['map/shape'][:])

        if has_names:
            names = loaded['names'][:]
        else:
            names = None

        return cmap, names


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
    # if not isinstance(data, np.ndarray):
    #     try:
    #         data = data.todense()
    #     except AttributeError:
    #         raise RuntimeError('Supplied data object does not appear to be a sparse matrix, nor a numpy array')

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
        if not isinstance(data, scipy.sparse.csr.csr_matrix):
            data = scipy.sparse.csr_matrix(data)
        save_h5map(fname, data, names)

    else:
        raise RuntimeError('Unimplemented format [{}]'.format(fmt))


def read_map(fname, fmt, delim=',', has_names=False, make_dense=False):
    """
    Read a contact map from file. The map is returned as a 2D
    numpy array. File formats are 'CSV' or a simple HDF5 with
    a single dataset named 'map'. In the case of a named CSV,
    a tuple will be returned (matrix, rownames, colnames).
    :param fname: input file name
    :param fmt: file format 'csv' or 'h5'
    :param delim: optional delimiter for csv
    :param has_names: for CSV, the first row and column specify names.
    :return: 
    """
    data = None
    cols = None

    if fmt == 'csv':
        with io_utils.open_input(fname) as csv_input:
            try:
                if has_names:
                    df = pandas.read_csv(csv_input, delimiter=delim, index_col=0, header=0)
                    if (df.index != df.columns).all():
                        raise IOError('row and column names do not agree')
                    cols = np.asarray(df.columns.tolist())
                else:
                    df = pandas.read_csv(csv_input, delimiter=delim, header=None)

            except ValueError as er:
                print 'Invalid data. Expected numeric values. Perhaps wrong delimiter used'
                raise er

            # always return floats
            data = df.as_matrix().astype(np.float)
            # pandas produces fortran arrays, we assume c-style
            data = np.ascontiguousarray(data)

    elif fmt == 'h5':
        data, cols = load_h5map(fname, has_names)
        if make_dense:
            data = np.ascontiguousarray(data.todense())

    else:
        raise RuntimeError('Unimplemented format [{}]'.format(fmt))

    return data, cols

