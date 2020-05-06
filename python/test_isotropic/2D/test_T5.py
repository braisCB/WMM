import sys
sys.path.append('../../')
import numpy as np
from compute_results import compute_results, aug_compute_results
from matplotlib import pyplot as plt


if __name__ == '__main__':

    h = 1
    tam = 101
    kappa = 1.
    tau = 8.

    initials = tam//2 * np.ones(2, dtype=int)

    gradient = np.empty((tam, tam, 2))
    ground_truth = np.empty((tam, tam))
    gradient_abs = np.empty((tam, tam))

    for i in range(tam):
        for j in range(tam):
            pos = h * (np.array([i, j]) - initials)
            zeta = np.sqrt(pos.dot(pos.T))
            ground_truth[i,j] = kappa * zeta - tau * np.sin(zeta / tau)
            gradient[i, j] = pos / zeta * (kappa - np.cos(zeta / tau))
            gradient_abs[i, j] = (kappa - np.cos(zeta / tau))

    gradient[np.isnan(gradient)] = 0.
    gradient_abs[np.isnan(gradient_abs)] = 0.

    x = np.arange(tam) - initials[1]
    y = np.arange(tam) - initials[0]
    xx, yy = np.meshgrid(x, y)
    plt.figure()
    plt.contour(xx, yy, ground_truth, 30)
    plt.title('$T_5$', fontsize=32)
    plt.show()

    inner_nodes_list = [3, 5]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
    for inner_nodes in inner_nodes_list:
        print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
