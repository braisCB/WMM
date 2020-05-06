import sys
sys.path.append('../../../')
import numpy as np
from compute_results import aug_compute_results


if __name__ == '__main__':

    h = 1.
    tam = 101

    initials = tam//2 * np.ones(2, dtype=int)

    gradient = np.empty((tam, tam, 2))
    ground_truth = np.empty((tam, tam))
    gradient_abs = np.empty((tam, tam))

    x = np.arange(tam) - initials[1]
    y = np.arange(tam) - initials[0]
    xx, yy = np.meshgrid(x, y)

    for i in range(tam):
        for j in range(tam):
            pos = h * (initials - np.array([i, j]))
            ground_truth[i,j] = np.square(pos[0]) / 25 + np.square(pos[1]) / 9
            gradient[i, j] = 2 * pos / np.array([25, 9])
            gradient_abs[i, j] = np.linalg.norm(gradient[i, j])

    inner_nodes_list = [3, 5]

    niter = 1

    for inner_nodes in inner_nodes_list:
        print('Computing MM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
