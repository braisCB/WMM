import numpy as np
from compute_results import compute_results, aug_compute_results


if __name__ == '__main__':

    h = [.1, .2, .1]
    tam = 41

    initials = tam//2 * np.ones(3, dtype=int)

    gradient = np.empty((tam, tam, tam, 3))
    ground_truth = np.empty((tam, tam, tam))
    gradient_abs = np.empty((tam, tam, tam))

    for i in range(tam):
        for j in range(tam):
            for k in range(tam):
                pos = h * (np.array([i, j, k]) - initials)
                norm = np.sqrt(pos.dot(pos.T))
                ground_truth[i,j,k] = norm
                gradient[i, j, k] = pos / norm if norm > 0 else [1, 0, 0]
                gradient_abs[i, j, k] = np.linalg.norm(gradient[i, j, k])

    # gradient[np.isnan(gradient)] = 0.
    # gradient_abs[np.isnan(gradient_abs)] = 0.

    inner_nodes_list = [1, 2]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
    for inner_nodes in inner_nodes_list:
        print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
