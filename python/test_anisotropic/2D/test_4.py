import sys
sys.path.append('../../')
import numpy as np
from compute_results_fp import compute_results, aug_compute_results
import mm
from matplotlib import pyplot as plt


if __name__ == '__main__':

    tam_gt = 385
    gradient_gt = np.zeros((tam_gt, tam_gt, 3))
    initials_gt = tam_gt // 2 * np.ones((1, 2), dtype=int)
    h_gt = 1 / (tam_gt - 1)

    A = 0.1225
    m = 2
    a = 0.5
    beta = 0

    for y in range(tam_gt):
        for x in range(tam_gt):
            p = (h_gt * (np.array([y, x]) - initials_gt))[0]
            y_i = A * np.sin(m * np.pi * p[1] / a + beta)
            grad_c = A * m * np.pi / a * np.cos(m * np.pi * p[1] / a + beta)
            if p[0] < -.25 + y_i:
                F1 = .2
                F2 = .8
            elif p[0] < 0. + y_i:
                F1 = 1.
                F2 = 1.
            elif p[0] < .25 + y_i:
                F1 = 1.
                F2 = 3.
            else:
                F1 = .2
                F2 = .8
            factor = np.sqrt(np.square(F2 / F1) - 1.) / np.sqrt(1. + np.square(grad_c))
            gradient_gt[y, x, 2] = F2
            gradient_gt[y, x, 0] = factor * grad_c
            gradient_gt[y, x, 1] = factor * -1
    ground_truth_gt = mm.aug_ani_wmm2D_fp(gradient_gt, initials_gt, h_gt, 15, 'linear', 'golden_search')
    # ground_truth_gt = mm.oum2D(gradient_gt, initials_gt, h_gt, 15, 'golden_search')
    plt.contour(ground_truth_gt, 50)
    plt.show()

    reductions = [2, 4, 8, 16]
    for reduction in reductions:
        tam = (tam_gt - 1) / reduction + 1
        initials = tam // 2 * np.ones((1, 2), dtype=int)
        h = reduction * h_gt

        gradient = gradient_gt[::reduction, ::reduction, :]

        ground_truth = ground_truth_gt[::reduction, ::reduction]

        inner_nodes_list = [3, 5, 15]

        niter = 1
        print('SIZE : ', tam)
        print('Computing WMM')
        compute_results(gradient=gradient, ground_truth=ground_truth,
                        initial_points=initials, niter=niter, h=h)
        for inner_nodes in inner_nodes_list:
            print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
            aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                                initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
