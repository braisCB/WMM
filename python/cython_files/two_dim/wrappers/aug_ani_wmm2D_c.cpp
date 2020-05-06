#include "mm2D_c.h"
#include "../../../../src/2D/aug_ani_wave_mm_2D.cpp"
#include "../../../../src/TYPES/WMMStructs.h"
#include <stdio.h>


void aug_ani_wmm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    double *h, int interp, int search, int N, int gamma, double *u_surface
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

    wmm::WmmAugAniSurface2D(image_grid, node_initials, hs, interp, search, N, gamma, u_surface_grid);

}

