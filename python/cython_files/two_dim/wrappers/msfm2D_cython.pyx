import numpy as np
cimport numpy as np

cdef extern from "mm2D_c.h":
    void msfm2D_c(
        double * image, int ncols, int nrows, int ndims, int * initials, int ninitials,
        double * h, int order, double * u_surface
    )

def msfm2d_cython(
        np.ndarray[np.double_t,ndim=3] image, np.ndarray[np.int32_t,ndim=2] initials,
        np.ndarray[np.double_t,ndim=1] h, int order
):
    image = np.ascontiguousarray(image)
    initials = np.ascontiguousarray(initials)
    h = np.ascontiguousarray(h)
    cdef int nrows = image.shape[0]
    cdef int ncols = image.shape[1]
    cdef int ndims = image.shape[2]
    cdef int ninitials = initials.shape[0]
    cdef np.ndarray[np.double_t,ndim=2,mode='c'] u_surface = np.Inf * np.ones((nrows, ncols))
    msfm2D_c(<double *> image.data, nrows, ncols, ndims, <int *> initials.data, ninitials,
            <double *> h.data, order, <double *> u_surface.data)
    return u_surface
