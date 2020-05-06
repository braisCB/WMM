import sys
sys.path.append('../../')
import numpy as np
import matplotlib.pyplot as plt
from compute_results import compute_results, aug_compute_results


if __name__ == '__main__':

    h = .5
    tam = 101
    degrees = [0, 21, 39, 57, 75]

    initials = tam//2 * np.ones(2, dtype=int)

    for degree in degrees:
        radian = degree * np.pi / 180.
        rotation_matrix = np.array([[np.cos(radian), np.sin(radian)], [-np.sin(radian), np.cos(radian)]])

        gradient = np.empty((tam, tam, 2))
        ground_truth = np.empty((tam, tam))
        gradient_abs = np.empty((tam, tam))

        for i in range(tam):
            for j in range(tam):
                pos = h * (initials - np.array([i, j]))
                pos = rotation_matrix @ pos
                ground_truth[i,j] = np.square(pos[0]) / 25 + np.square(pos[1]) / 9
                gradient[i, j] = rotation_matrix.T @ (2 * pos / np.array([25, 9]))
                gradient_abs[i, j] = np.linalg.norm(gradient[i, j])

        #plt.quiver(gradient[:, :, 1], gradient[:, :, 0], angles='xy')
        #plt.show()

        inner_nodes_list = [3, 5]

        niter = 1

        print('Computing WMM')
        compute_results(gradient=gradient, ground_truth=ground_truth,
                        initial_points=initials, niter=niter, h=h)
        for inner_nodes in inner_nodes_list:
            print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
            aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                                initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
