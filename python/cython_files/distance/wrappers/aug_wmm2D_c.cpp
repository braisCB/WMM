#include "distance2D_c.h"
#include "../../../../src/TYPES/DistanceWMMStructs.h"
#include "../../../../src/DistanceSurface/aug_wave_mm_2D.cpp"
#include <stdio.h>


void aug_wmm2D_c(
    double *image, int ncols, int nrows, int ndims, int *initials, int ninitials,
    int *finals, int nfinals,
    double *h, int interp, int search, int N, int gamma, double *u_surface
) {

    wmm::Grid image_grid(image, nrows, ncols, ndims, 'C');
    std::vector<wmm::Node> node_initials, node_finals;
    wmm::Node p;
    for (int i = 0; i < ninitials; i++) {
        p = wmm::Node(initials[2*i], initials[2*i + 1]);
        node_initials.push_back(p);
    }
    for (int i = 0; i < nfinals; i++) {
        p = wmm::Node(finals[2*i], finals[2*i + 1]);
        node_finals.push_back(p);
    }
    wmm::NodeD hs = wmm::NodeD(h[0], h[1]);
    wmm::Grid u_surface_grid(u_surface, nrows, ncols, 2, 'C');

    wmm::WmmAugSurface2D(image_grid, node_initials, node_finals, hs, interp, search, N, gamma, u_surface_grid);

}

