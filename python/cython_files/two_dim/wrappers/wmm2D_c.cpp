#include "mm2D_c.h"
#include "../../../../src/2D/wave_mm_2D.cpp"
#include "../../../../src/TYPES/WMMStructs.h"
#include <stdio.h>


void wmm2D_c(
    double *image, int nrows, int ncols, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, double *u_surface
) {

    wmm::Grid image_grid(image, nrows, ncols, ndims, 'C');
    std::vector<wmm::Node> node_initials;
    wmm::Node p;
    for (int i = 0; i < ninitials; i++) {
        p = wmm::Node(initials[2*i], initials[2*i + 1]);
        node_initials.push_back(p);
    }
    wmm::NodeD hs = wmm::NodeD(h[0], h[1]);
    wmm::Grid u_surface_grid(u_surface, nrows, ncols, 1, 'C');

    switch (interp) {
        case wmm::I_LINEAR:
            wmm::WmmIsoSurface2D<1>(image_grid, node_initials, hs, interp, search, u_surface_grid);
            break;
        case wmm::I_QUADATRIC:
            wmm::WmmIsoSurface2D<5>(image_grid, node_initials, hs, interp, search, u_surface_grid);
            break;
        default: //I_SPLINE, I_HERMITE and I_PCHIP
            wmm::WmmIsoSurface2D<4>(image_grid, node_initials, hs, interp, search, u_surface_grid);
    }

}

