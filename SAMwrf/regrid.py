"""
Sparse matrix multiplication (SMM) using scipy.sparse library.
"""

import xarray as xr
import scipy.sparse as sps
import warnings


def read_weights(filename, n_in, n_out):
    '''
    Read regridding weights into a scipy sparse COO matrix.
    Parameters
    ----------
    filename : str
        Offline weight file generated by ESMPy.
    N_in, N_out : integers
        ``(N_out, N_in)`` will be the shape of the returning sparse matrix.
        They are the total number of grid boxes in input and output grids::
              N_in = Nx_in * Ny_in
              N_out = Nx_out * Ny_out
        We need them because the shape cannot always be infered from the
        largest column and row indices, due to unmapped grid boxes.
    Returns
    -------
    A : scipy sparse COO matrix.
    '''
    ds_w = xr.open_dataset(filename)

    col = ds_w['col'].values - 1  # Python starts with 0
    row = ds_w['row'].values - 1
    S = ds_w['S'].values

    weights = sps.coo_matrix((S, (row, col)), shape=[n_out, n_in])
    return weights


def apply_weights(weights, indata, shape_in, shape_out):
    '''
    Apply regridding weights to data.
    Parameters
    ----------
    A : scipy sparse COO matrix
    indata : numpy array of shape ``(..., n_lat, n_lon)`` or ``(..., n_y, n_x)``.
        Should be C-ordered. Will be then tranposed to F-ordered.
    shape_in, shape_out : tuple of two integers
        Input/output data shape for unflatten operation.
        For rectilinear grid, it is just ``(n_lat, n_lon)``.
    Returns
    -------
    outdata : numpy array of shape ``(..., shape_out[0], shape_out[1])``.
        Extra dimensions are the same as `indata`.
        If input data is C-ordered, output will also be C-ordered.
    '''

    # COO matrix is fast with F-ordered array but slow with C-array, so we
    # take in a C-ordered and then transpose)
    # (CSR or CRS matrix is fast with C-ordered array but slow with F-array)
    if not indata.flags['C_CONTIGUOUS']:
        warnings.warn("Input array is not C_CONTIGUOUS. "
                      "Will affect performance.")

    # get input shape information
    shape_horiz = indata.shape[-2:]
    extra_shape = indata.shape[0:-2]

    assert shape_horiz == shape_in, (
        'The horizontal shape of input data is {}, different from that of'
        'the regridder {}!'.format(shape_horiz, shape_in)
        )

    assert shape_in[0] * shape_in[1] == weights.shape[1], (
        "ny_in * nx_in should equal to weights.shape[1]")

    assert shape_out[0] * shape_out[1] == weights.shape[0], (
        "ny_out * nx_out should equal to weights.shape[0]")

    # use flattened array for dot operation
    indata_flat = indata.reshape(-1, shape_in[0]*shape_in[1])
    outdata_flat = weights.dot(indata_flat.T).T

    # unflattened output array
    outdata = outdata_flat.reshape(
        [*extra_shape, shape_out[0], shape_out[1]])
    return outdata
