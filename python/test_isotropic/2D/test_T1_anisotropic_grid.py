import sys
sys.path.append('../../')
import numpy as np
from compute_results import compute_results, aug_compute_results


if __name__ == '__main__':

    h = np.array([.2, .1])
    tam = 101

    initials = tam//2 * np.ones((1,2), dtype=int)

    gradient = np.empty((tam, tam, 2))
    ground_truth = np.empty((tam, tam))
    gradient_abs = np.empty((tam, tam))

    for i in range(tam):
        for j in range(tam):
            diff = h * (np.array([i, j]) - initials[0])
            norm = np.sqrt(diff.dot(diff.T))
            ground_truth[i,j] = norm
            if norm == 0:
                diff = np.array([1., 0.])
                norm = 1
            gradient[i, j] = diff / norm
            gradient_abs[i, j] = np.linalg.norm(gradient[i, j])

    inner_nodes_list = [3, 5]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
    for inner_nodes in inner_nodes_list:
        print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
