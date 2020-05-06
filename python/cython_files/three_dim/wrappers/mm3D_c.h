// wmm2D.h: numpy arrays from cython , double*

void wmm3D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, double *u_surface
);

void aug_wmm3D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, int N, int gamma, double *u_surface
);

void ani_wmm3D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, double *u_surface
);

void ani_wmm3D_riemann_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, double *u_surface, int *dir_surface
);

void aug_ani_wmm3D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, int N, int gamma, double *u_surface
);
