import sys
sys.path.append('../../')
import numpy as np
from compute_results import compute_results, aug_compute_results, compute_oum_results
import mm
from matplotlib import pyplot as plt


if __name__ == '__main__':

    tam_gt = 385
    gradient_gt = np.zeros((tam_gt, tam_gt, 2))
    initials_gt = tam_gt // 2 * np.ones((1, 2), dtype=int)
    h_gt = 1 / (tam_gt - 1)

    for y in range(tam_gt):
        for x in range(tam_gt):
            p = (h_gt * (np.array([y, x]) - initials_gt))[0]
            gradient_gt[y, x, 0] = 0.9 * 2 * np.pi * np.cos(2 * np.pi * p[0]) * np.sin(2 * np.pi * p[1])
            gradient_gt[y, x, 1] = 0.9 * 2 * np.pi * np.sin(2 * np.pi * p[0]) * np.cos(2 * np.pi * p[1])
    # ground_truth_gt = mm.ani_wmm2D(gradient_gt, initials_gt, h_gt, 'pchip', 'golden_search')
    ground_truth_gt = mm.oum2D(gradient_gt, initials_gt, h_gt, 6, 'golden_search')
    plt.contour(ground_truth_gt, 30)
    plt.show()

    reductions = [2, 4, 8, 16]
    max_reduction = np.max(reductions)
    for reduction in reductions:
        tam = (tam_gt - 1) / reduction + 1
        initials = tam // 2 * np.ones((1, 2), dtype=int)
        h = reduction * h_gt

        gradient = gradient_gt[::reduction, ::reduction, :]

        ground_truth = ground_truth_gt[::reduction, ::reduction]

        inner_nodes_list = [3, 5]

        niter = 10
        print('SIZE : ', tam)
        print('Computing WMM')
        compute_oum_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, h=h, gamma=6)
        compute_results(gradient=gradient, ground_truth=ground_truth,
                        initial_points=initials, niter=niter, h=h)
        for inner_nodes in inner_nodes_list:
            print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
            aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                                initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
