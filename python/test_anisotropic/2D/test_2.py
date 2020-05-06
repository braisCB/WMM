import sys
sys.path.append('../../')
import numpy as np
from compute_results import compute_results, aug_compute_results
import mm
from matplotlib import pyplot as plt


if __name__ == '__main__':

    h = 1
    tam = 101
    extra_dim = 5

    initials = tam//2 * np.ones((1,2), dtype=int)

    gradient = np.zeros((tam, tam, 2))
    gradient[:,:,1] = np.sqrt(2.)

    tam_gt = (tam - 1) * extra_dim + 1
    gradient_gt = np.zeros((tam_gt, tam_gt, 2))
    gradient_gt[:,:,1] = np.sqrt(2.)
    initials_gt = tam_gt//2 * np.ones((1,2), dtype=int)
    h_gt = h / extra_dim

    ground_truth = mm.ani_wmm2D(gradient_gt, initials_gt, h_gt, 'pchip', 'golden_search')
    plt.contour(ground_truth)
    plt.show()
    ground_truth = ground_truth[::extra_dim, ::extra_dim]

    inner_nodes_list = [3, 5]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
    for inner_nodes in inner_nodes_list:
        print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
