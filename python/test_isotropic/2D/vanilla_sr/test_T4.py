import sys
sys.path.append('../../../')
import numpy as np
from compute_results import aug_compute_results
from matplotlib import pyplot as plt


if __name__ == '__main__':

    h = 1
    tam = 101
    kappa = 1.
    tau = 2.

    initials = tam//2 * np.ones(2, dtype=int)

    gradient = np.zeros((tam, tam, 2))
    ground_truth = np.empty((tam, tam))
    gradient_abs = np.empty((tam, tam))

    for i in range(tam):
        for j in range(tam):
            pos = h * (initials - np.array([i, j]))
            zeta = np.sqrt(pos.dot(pos.T))
            ground_truth[i, j] = kappa * zeta - tau * np.sin(zeta / tau)
            gradient[i, j] = pos / zeta * (kappa - np.cos(zeta / tau))
            gradient_abs[i, j] = (kappa - np.cos(zeta / tau))

    gradient[np.isnan(gradient)] = 0.
    gradient_abs[np.isnan(gradient_abs)] = 0.

    inner_nodes_list = [3, 5]

    niter = 1

    for inner_nodes in inner_nodes_list:
        print('Computing MM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)

