// wmm2D.h: numpy arrays from cython , double*

void wmm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    int *finals, int nfinals,
    double *h, int interp, int search, double *u_surface
);

void aug_wmm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    int *finals, int nfinals,
    double *h, int interp, int search, int N, int gamma, double *u_surface
);

void msfm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    int *finals, int nfinals,
    double *h, int order, double *u_surface
);

void fmm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    int *finals, int nfinals,
    double *h, int order, double *u_surface
);

