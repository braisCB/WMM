import numpy as np
from cython_files.two_dim.wmm2D_cython import wmm2d_cython
from cython_files.two_dim.wmm2D_cartesian_cython import wmm2d_cartesian_cython
from cython_files.two_dim.ani_wmm2D_cython import ani_wmm2d_cython
from cython_files.two_dim.ani_wmm2D_fp_cython import ani_wmm2D_fp_cython
from cython_files.two_dim.aug_wmm2D_cython import aug_wmm2d_cython
from cython_files.two_dim.aug_ani_wmm2D_cython import aug_ani_wmm2d_cython
from cython_files.two_dim.aug_ani_wmm2D_fp_cython import aug_ani_wmm2D_fp_cython
from cython_files.two_dim.aug_ani_wmm2D_riemann_cython import aug_ani_wmm2D_riemann_cython
from cython_files.two_dim.ani_wmm2D_riemann_cython import ani_wmm2D_riemann_cython
from cython_files.two_dim.oum2D_cython import oum2d_cython
from cython_files.two_dim.fmm2D_cython import fmm2d_cython
from cython_files.two_dim.msfm2D_cython import msfm2d_cython
from cython_files.three_dim.wmm3D_cython import wmm3d_cython
from cython_files.three_dim.aug_wmm3D_cython import aug_wmm3d_cython
from cython_files.three_dim.ani_wmm3D_cython import ani_wmm3d_cython
#from cython_files.three_dim.ani_wmm3D_riemann_cython import ani_wmm3D_riemann_cython
from cython_files.three_dim.aug_ani_wmm3D_cython import aug_ani_wmm3d_cython


def wmm2D(image, initials, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (0, 0)], 'edge')
    image = np.asarray(image, dtype=np.double)

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
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

    u_surface = wmm2d_cython(image, initials, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad]


def wmm2D_cartesian(image, initials, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (0, 0)], 'edge')
    image = np.asarray(image, dtype=np.double)

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
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

    u_surface = wmm2d_cython(image, initials, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad]


def ani_wmm2D(image, initials, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (0, 0)], 'edge')
    image = np.asarray(image, dtype=np.double)

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
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

    u_surface = ani_wmm2d_cython(image, initials, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad]


def ani_wmm2D_fp(image, initials, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (0, 0)], 'edge')
    image = np.asarray(image, dtype=np.double)

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
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

    u_surface = ani_wmm2D_fp_cython(image, initials, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad]


def ani_wmm2D_riemann(image, initials, h, interp, search, pad=0):

    image = np.asarray(image, dtype=np.double)

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)

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

    u_surface = ani_wmm2D_riemann_cython(image, initials, h, interp_int, search_int)

    return u_surface


def aug_wmm2D(image, initials, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    u_surface = aug_wmm2d_cython(image, initials, h, interp_int, search_int, N, gamma)

    return u_surface


def aug_ani_wmm2D(image, initials, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    u_surface = aug_ani_wmm2d_cython(image, initials, h, interp_int, search_int, N, gamma)

    return u_surface


def aug_ani_wmm2D_fp(image, initials, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    u_surface = aug_ani_wmm2D_fp_cython(image, initials, h, interp_int, search_int, N, gamma)

    return u_surface


def aug_ani_wmm2D_riemann(image, initials, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    u_surface = aug_ani_wmm2D_riemann_cython(image, initials, h, interp_int, search_int, N, gamma)

    return u_surface


def oum2D(image, initials, h, radius, search):

    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    if search == 'gradient':
        search_int = 0
    elif search == 'hopf_lax':
        search_int = 1
    elif search == 'golden_search':
        search_int = 2
    else:
        raise Exception('incorrect search value')

    u_surface = oum2d_cython(image, initials, h, radius, search_int)

    return u_surface


def fmm2D(image, initials, h, order):

    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    u_surface = fmm2d_cython(image, initials, h, order)

    return u_surface


def msfm2D(image, initials, h, order):

    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
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

    u_surface = msfm2d_cython(image, initials, h, order)

    return u_surface


def wmm3D(image, initials, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (pad, pad), (0,0)], 'edge')
    image = np.asarray(image, dtype=np.double)
    # augmented_image = np.augmented_image(image, dtype=np.double)
    # image = np.ascontiguousarray(image)
    # image = image.copy()

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 2] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1]) &
        (initials[:, 2] >= image.shape[2])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 3

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 3
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if interp == 'bilinear':
        interp_int = 0
    elif interp == 'quad':
        interp_int = 4
    elif interp == 'pchip':
        interp_int = 2
    elif interp == 'hermite':
        interp_int = 1
    elif interp == 'spline':
        interp_int = 3
    elif interp == 'cubic':
        interp_int = 9
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

    u_surface = wmm3d_cython(image, initials, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad, pad:-pad]


def ani_wmm3D(image, initials, h, interp, search, pad=3):

    image = np.pad(image, [(pad, pad), (pad, pad), (pad, pad), (0,0)], 'edge')
    image = np.asarray(image, dtype=np.double)
    # augmented_image = np.augmented_image(image, dtype=np.double)
    # image = np.ascontiguousarray(image)
    # image = image.copy()

    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    initials += pad
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 2] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1]) &
        (initials[:, 2] >= image.shape[2])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 3

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 3
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if interp == 'bilinear':
        interp_int = 0
    elif interp == 'quad':
        interp_int = 4
    elif interp == 'pchip':
        interp_int = 2
    elif interp == 'hermite':
        interp_int = 1
    elif interp == 'spline':
        interp_int = 3
    elif interp == 'cubic':
        interp_int = 9
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

    u_surface = ani_wmm3d_cython(image, initials, h, interp_int, search_int)

    return u_surface[pad:-pad, pad:-pad, pad:-pad]


# def ani_wmm3D_riemann(image, initials, h, interp, search):
#
#     # image = np.pad(image, [(pad, pad), (pad, pad), (pad, pad), (0,0)], 'edge')
#     image = np.asarray(image, dtype=np.double)
#     # augmented_image = np.augmented_image(image, dtype=np.double)
#     # image = np.ascontiguousarray(image)
#     # image = image.copy()
#
#     initials = np.asarray(initials, dtype=np.int32)
#     if initials.ndim == 1:
#         initials = np.expand_dims(initials, axis=0)
#     # initials += pad
#     if isinstance(h, float) or isinstance(h, int):
#         h = [h, h, h]
#     h = np.asarray(h, dtype=np.double)
#
#     check_initials = len(np.where(
#         (initials[:, 0] < 0) &
#         (initials[:, 1] < 0) &
#         (initials[:, 2] < 0) &
#         (initials[:, 0] >= image.shape[0]) &
#         (initials[:, 1] >= image.shape[1]) &
#         (initials[:, 2] >= image.shape[2])
#     )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 3
#
#     if not check_initials:
#         raise Exception('initial points do not lie inside the image')
#
#     check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 3
#     if not check_h:
#         raise Exception('h has an incorrect value or shape')
#
#     if interp == 'bilinear':
#         interp_int = 0
#     elif interp == 'quad':
#         interp_int = 4
#     elif interp == 'pchip':
#         interp_int = 2
#     elif interp == 'hermite':
#         interp_int = 1
#     elif interp == 'spline':
#         interp_int = 3
#     elif interp == 'cubic':
#         interp_int = 9
#     else:
#         raise Exception('incorrect interp value')
#
#     if search == 'gradient':
#         search_int = 0
#     elif search == 'hopf_lax':
#         search_int = 1
#     elif search == 'golden_search':
#         search_int = 2
#     else:
#         raise Exception('incorrect search value')
#
#     u_surface, propagation = ani_wmm3D_riemann_cython(image, initials, h, interp_int, search_int)
#
#     return u_surface, propagation # [pad:-pad, pad:-pad, pad:-pad]


def aug_wmm3D(image, initials, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 2] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1]) &
        (initials[:, 2] >= image.shape[2])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 3

    if not isinstance(gamma, int):
        raise Exception('gamma should be an integer')

    if gamma < 0:
        raise Exception('gamma should be positive')

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 3
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if interp == 'bilinear':
        interp_int = 0
        N = 1
    elif interp == 'quad':
        interp_int = 4
        N = 4
    elif interp == 'pchip':
        interp_int = 2
        N = 11
    elif interp == 'hermite':
        interp_int = 1
        N = 11
    elif interp == 'spline':
        interp_int = 3
        N = 11
    elif interp == 'cubic':
        interp_int = 9
        N = 9
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

    u_surface = aug_wmm3d_cython(image, initials, h, interp_int, search_int, N, gamma)

    return u_surface


def aug_ani_wmm3D(image, initials, h, gamma, interp, search):

    gamma += 1
    image = np.asarray(image, dtype=np.double)
    initials = np.asarray(initials, dtype=np.int32)
    if initials.ndim == 1:
        initials = np.expand_dims(initials, axis=0)
    if isinstance(h, float) or isinstance(h, int):
        h = [h, h, h]
    h = np.asarray(h, dtype=np.double)

    check_initials = len(np.where(
        (initials[:, 0] < 0) &
        (initials[:, 1] < 0) &
        (initials[:, 2] < 0) &
        (initials[:, 0] >= image.shape[0]) &
        (initials[:, 1] >= image.shape[1]) &
        (initials[:, 2] >= image.shape[2])
    )[0]) == 0 and initials.ndim == 2 and initials.shape[1] == 3

    if not isinstance(gamma, int):
        raise Exception('gamma should be an integer')

    if gamma < 0:
        raise Exception('gamma should be positive')

    if not check_initials:
        raise Exception('initial points do not lie inside the image')

    check_h = len(np.where(h <= 0)[0]) == 0 and h.ndim == 1 and h.shape[0] == 3
    if not check_h:
        raise Exception('h has an incorrect value or shape')

    if interp == 'bilinear':
        interp_int = 0
        N = 1
    elif interp == 'quad':
        interp_int = 4
        N = 4
    elif interp == 'pchip':
        interp_int = 2
        N = 11
    elif interp == 'hermite':
        interp_int = 1
        N = 11
    elif interp == 'spline':
        interp_int = 3
        N = 11
    elif interp == 'cubic':
        interp_int = 9
        N = 9
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

    u_surface = aug_ani_wmm3d_cython(image, initials, h, interp_int, search_int, N, gamma)

    return u_surface

