import sys
sys.path.append('../../../')
import numpy as np
from compute_results import aug_compute_results


if __name__ == '__main__':

    h = 1
    tam = 101

    initials = tam//2 * np.ones(2, dtype=int)

    gradient = np.empty((tam, tam, 2))
    ground_truth = np.empty((tam, tam))
    gradient_abs = np.empty((tam, tam))

    for i in range(tam):
        for j in range(tam):
            pos = h * (np.array([i, j]) - initials)
            ground_truth[i,j] = np.square(pos[0]) / 100 + np.square(pos[1]) / 20
            gradient[i, j] = 2 * pos / np.array([100, 20])
            gradient_abs[i, j] = np.linalg.norm(gradient[i, j])

    # gradient[np.isnan(gradient)] = 0.
    # gradient_abs[np.isnan(gradient_abs)] = 0.

    inner_nodes_list = [3, 5]

    niter = 1

    for inner_nodes in inner_nodes_list:
        print('Computing MM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
