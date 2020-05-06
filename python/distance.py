import numpy as np
from cython_files.distance.wmm2D_cython import wmm2d_cython
from cython_files.distance.aug_wmm2D_cython import aug_wmm2d_cython
from cython_files.distance.fmm2D_cython import fmm2d_cython
from cython_files.distance.msfm2D_cython import msfm2d_cython


def wmm2D(image, initials, finals, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (0, 0)], 'edge')
    image = np.asarray(image, dtype=np.double)

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
    finals = np.asarray(finals, dtype=np.int32)
    if finals.ndim == 1:
        finals = np.expand_dims(finals, axis=0)
    finals += pad
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 2

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 2
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if interp == 'linear':
        interp_int = 0
    elif interp == 'quad':
        interp_int = 4
    elif interp == 'pchip':
        interp_int = 2
    elif interp == 'hermite':
        interp_int = 1
    elif interp == 'spline':
        interp_int = 3
    else:
        raise Exception('incorrect interp value')

    if search == 'gradient':
        search_int = 0
    elif search == 'hopf_lax':
        search_int = 1
    elif search == 'golden_search':
        search_int = 2
    else:
        raise Exception('incorrect search value')

    u_surface = wmm2d_cython(image, initials, finals, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad, 0], u_surface[pad:-pad, pad:-pad, 1]


def aug_wmm2D(image, initials, finals, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    finals = np.asarray(finals, dtype=np.int32)
    if finals.ndim == 1:
        finals = np.expand_dims(finals, axis=0)
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 2

    if not isinstance(gamma, int):
        raise Exception('gamma should be an integer')

    if gamma <= 0:
        raise Exception('gamma should be positive')

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 2
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if interp == 'linear':
        interp_int = 0
        N = 0
    elif interp == 'quad':
        interp_int = 4
        N = 3
    elif interp in ['pchip', 'hermite']:
        interp_int = 2 if interp == 'pchip' else 1
        N = 2
    elif interp == 'spline':
        interp_int = 3
        N = 2
    else:
        raise Exception('incorrect interp value')

    if search == 'gradient':
        search_int = 0
    elif search == 'hopf_lax':
        search_int = 1
    elif search == 'golden_search':
        search_int = 2
    else:
        raise Exception('incorrect search value')

    u_surface = aug_wmm2d_cython(image, initials, finals, h, interp_int, search_int, N, gamma)

    return u_surface[..., 0], u_surface[..., 1]


def fmm2D(image, initials, finals, h, order):

    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    finals = np.asarray(finals, dtype=np.int32)
    if finals.ndim == 1:
        finals = np.expand_dims(finals, axis=0)
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 2

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 2
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if order not in [1, 2]:
        raise Exception('incorrect order value: should be 1 or 2')

    u_surface = fmm2d_cython(image, initials, finals, h, order)

    return u_surface[..., 0], u_surface[..., 1]


def msfm2D(image, initials, finals, h, order):

    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    finals = np.asarray(finals, dtype=np.int32)
    if finals.ndim == 1:
        finals = np.expand_dims(finals, axis=0)
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 2

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 2
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if order not in [1, 2]:
        raise Exception('incorrect order value: should be 1 or 2')

    u_surface = msfm2d_cython(image, initials, finals, h, order)

    return u_surface[..., 0], u_surface[..., 1]
