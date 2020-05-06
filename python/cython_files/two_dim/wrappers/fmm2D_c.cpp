#include "mm2D_c.h"
#include "../../../../src/2D/fmm_2D.cpp"
#include "../../../../src/TYPES/WMMStructs.h"
#include <stdio.h>


void fmm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int order, double *u_surface
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

    wmm::FMMIsoSurface2D(image_grid, node_initials, hs, order, u_surface_grid);

}

