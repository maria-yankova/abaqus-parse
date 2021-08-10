import numpy as np


def flatten_list(a):
    return [item for sublist in a for item in sublist]


def format_args_check(**kwargs):
    """
    Check types of parameters used in `format_arr`, 'format_list' and
    'format_dict' functions.

    """

    if 'depth' in kwargs and not isinstance(kwargs['depth'], int):
        raise ValueError('`depth` must be an integer.')

    if 'indent' in kwargs and not isinstance(kwargs['indent'], str):
        raise ValueError('`indent` must be a string.')

    if 'col_delim' in kwargs and not isinstance(kwargs['col_delim'], str):
        raise ValueError('`col_delim` must be a string.')

    if 'row_delim' in kwargs and not isinstance(kwargs['row_delim'], str):
        raise ValueError('`row_delim` must be a string.')

    if 'dim_delim' in kwargs and not isinstance(kwargs['dim_delim'], str):
        raise ValueError('`dim_delim` must be a string.')

    if 'format_spec' in kwargs and not isinstance(kwargs['format_spec'],
                                                  (str, list)):
        raise ValueError('`format_spec` must be a string or list of strings.')

    if 'assign' in kwargs:

        if not isinstance(kwargs['assign'], str):
            raise ValueError('`assign` must be a string.')


def format_arr(arr, depth=0, indent='\t', col_delim='\t', row_delim='\n',
               dim_delim='\n', format_spec='{}'):
    """
    Get a string representation of a Numpy array, formatted with indents.

    Parameters
    ----------
    arr : ndarray or list of ndarray
        Array of any shape to format as a string, or list of arrays whose
        shapes match except for the final dimension, in which case the arrays
        will be formatted horizontally next to each other.
    depth : int, optional
        The indent depth at which to begin the formatting.
    indent : str, optional
        The string used as the indent. The string which indents each line of
        the array is equal to (`indent` * `depth`).
    col_delim : str, optional
        String to delimit columns (the innermost dimension of the array).
        Default is tab character, \t.
    row_delim : str, optional
        String to delimit rows (the second-innermost dimension of the array).
        Default is newline character, \n.
    dim_delim : str, optional
        String to delimit outer dimensions. Default is newline character, \n.
    format_spec : str or list of str, optional
        Format specifier for the array or a list of format specifiers, one for 
        each array listed in `arr`.

    Returns
    -------
    str

    """

    # Validation:
    format_args_check(depth=depth, indent=indent, col_delim=col_delim,
                      row_delim=row_delim, dim_delim=dim_delim,
                      format_spec=format_spec)

    if isinstance(arr, np.ndarray):
        arr = [arr]

    out_shape = list(set([i.shape[:-1] for i in arr]))

    if len(out_shape) > 1:
        raise ValueError('Array shapes must be identical apart from the '
                         'innermost dimension.')

    if not isinstance(arr, (list, np.ndarray)):
        raise ValueError('Cannot format as array, object is '
                         'not an array or list of arrays: type is {}'.format(
                             type(arr)))

    if isinstance(format_spec, str):
        format_spec = [format_spec] * len(arr)

    elif isinstance(format_spec, list):

        fs_err_msg = ('`format_spec` must be a string or list of N strings '
                      'where N is the number of arrays specified in `arr`.')

        if not all([isinstance(i, str)
                    for i in format_spec]) or len(format_spec) != len(arr):
            raise ValueError(fs_err_msg)

    arr_list = arr
    out = ''
    dim_seps = ''
    d = arr_list[0].ndim

    if d == 1:
        out += (indent * depth)

        for sa_idx, sub_arr in enumerate(arr_list):
            for idx_col, col in enumerate(sub_arr):
                if idx_col == 0 and sa_idx == 0 and idx_col == len(sub_arr)-1 and sa_idx == len(arr_list)-1:
                    out += '   ' + format_spec[sa_idx].format(col)
                elif idx_col == len(sub_arr)-1 and sa_idx == len(arr_list)-1:
                    out += format_spec[sa_idx].format(col)
                elif idx_col == 0 and sa_idx == 0:
                    out += '   ' + format_spec[sa_idx].format(col) + col_delim
                else:
                    out += format_spec[sa_idx].format(col) + col_delim

        out += row_delim

    else:

        if d > 2:
            dim_seps = dim_delim * (d - 2)

        sub_arr = []
        for i in range(out_shape[0][0]):

            sub_arr_lst = []
            for j in arr_list:
                sub_arr_lst.append(j[i])

            sub_arr.append(format_arr(sub_arr_lst, depth, indent, col_delim,
                                      row_delim, dim_delim, format_spec))

        out = dim_seps.join(sub_arr)

    return out


def circle_points(r, n, centre):
    t = np.linspace(0, 2*np.pi, n)
    x = r * np.cos(t) + centre[0]
    y = r * np.sin(t) + centre[1]

    return np.c_[x, y]


def polar2cart_2D(r, θ):
    """
    Convert 2D polar coordinates to Cartesian coordinates.

    """

    x = r * np.cos(θ)
    y = r * np.sin(θ)

    return x, y


def order_coplanar_points(points, normal, anticlockwise=True):
    """
    Find the clockwise or anticlockwise ordering of a set of coplanar 3D points.

    Parameters
    ----------
    points : ndarray of shape (3, N)
        The set of coplanar points (three-vector columns) whose ordering is to be found.
    normal : ndarray of shape (3, 1)
        Column three-vector representing the normal of the plane on which all points lie.

    Returns
    -------
    Ordered indices of points according to a clockwise or anticlockwise direction when
    looking in the opposite direction to `normal`.

    """

    # Normalise `normal` to a unit vector:
    normal = normal / np.linalg.norm(normal)

    # Compute the centroid:
    centroid = np.mean(points, axis=1)

    # Get the direction vectors from each point to the centroid
    p = points - centroid[:, np.newaxis]

    # Use the first point as the reference point
    # Find cross product of each point with the reference point
    crs = np.cross(p[:, 0:1], p, axis=0)

    # Find the scalar triple product of the point pairs and the normal vector
    stp = np.einsum('ij,ik->k', normal, crs)

    # Find the dot product of the point pairs
    dot = np.einsum('ij,ik->k', p[:, 0:1], p)

    # Find signed angles from reference point to each other point
    ang = np.arctan2(stp, dot)
    ang_order = np.argsort(ang)

    if not anticlockwise:
        ang_order = ang_order[::-1]

    return ang_order


def search_keep_order(A, B):
    sort_idx = A.argsort()
    out = sort_idx[np.searchsorted(A, B, sorter = sort_idx)]
    idx_layer_0 = np.nonzero(B[:,None] == A)[1]
    
    return idx_layer_0