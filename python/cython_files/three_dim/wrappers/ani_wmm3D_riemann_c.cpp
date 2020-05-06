#include "mm3D_c.h"
#include "../../../../src/3D/ani_wave_mm_3D_riemann.cpp"
#include "../../../../src/TYPES/WMMStructs3D.h"
#include <stdio.h>


void ani_wmm3D_riemann_c(
    double *image, int nrows, int ncols, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, double *u_surface, int *dir_surface
) {

    wmm3D::Grid3 image_grid(image, nrows, ncols, ndims, 9, 'C');
    std::vector<wmm3D::Node3> node_initials;
    wmm3D::Node3 p;
    for (int i = 0; i < ninitials; i++) {
        p = wmm3D::Node3(initials[3*i], initials[3*i + 1], initials[3*i + 2]);
        node_initials.push_back(p);
    }
    wmm3D::Node3D hs = wmm3D::Node3D(h[0], h[1], h[2]);
    wmm3D::Grid3 u_surface_grid(u_surface, nrows, ncols, ndims, 1, 'C');
    wmm3D::Grid3_<int> dir_surface_grid(dir_surface, nrows, ncols, ndims, 3, 'C');

    switch (interp) {
        case wmm3D::I_BILINEAR:
            wmm3D::wmm3DAniSurface3D<0>(image_grid, node_initials, hs, interp, search, u_surface_grid, dir_surface_grid);
            break;
        case wmm3D::I_QUADATRIC:
            wmm3D::wmm3DAniSurface3D<4>(image_grid, node_initials, hs, interp, search, u_surface_grid, dir_surface_grid);
            break;
        case wmm3D::I_CUBIC:
            wmm3D::wmm3DAniSurface3D<9>(image_grid, node_initials, hs, interp, search, u_surface_grid, dir_surface_grid);
            break;
        default: //I_SPLINE, I_HERMITE and I_PCHIP
            wmm3D::wmm3DAniSurface3D<11>(image_grid, node_initials, hs, interp, search, u_surface_grid, dir_surface_grid);
    }

}

