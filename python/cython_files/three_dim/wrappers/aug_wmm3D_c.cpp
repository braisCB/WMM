#include "mm3D_c.h"
#include "../../../../src/3D/aug_wave_mm_3D.cpp"
#include "../../../../src/TYPES/WMMStructs3D.h"
#include <stdio.h>


void aug_wmm3D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, int N, int gamma, double *u_surface
) {

    wmm3D::Grid3 image_grid(image, nrows, ncols, ndims, 3, 'C');
    std::vector<wmm3D::Node3> node_initials;
    wmm3D::Node3 p;
    for (int i = 0; i < ninitials; i++) {
        p = wmm3D::Node3(initials[3*i], initials[3*i + 1], initials[3*i + 2]);
        node_initials.push_back(p);
    }
    wmm3D::Node3D hs = wmm3D::Node3D(h[0], h[1], h[2]);
    wmm3D::Grid3 u_surface_grid(u_surface, nrows, ncols, ndims, 1, 'C');

    wmm3D::augWmm3DIsoSurface3D(image_grid, node_initials, hs, interp, search, N, gamma, u_surface_grid);

}

