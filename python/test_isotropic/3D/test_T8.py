import numpy as np
from compute_results import compute_results, aug_compute_results
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


if __name__ == '__main__':

    h = 1
    tam = 51

    initials = tam//2 * np.ones(3, dtype=int)

    gradient = np.empty((tam, tam, tam, 3))
    ground_truth = np.empty((tam, tam, tam))
    gradient_abs = np.empty((tam, tam, tam))

    factor = np.array([100, 20, 20])

    for i in range(tam):
        for j in range(tam):
            for k in range(tam):
                pos = h * (np.array([i, j, k]) - initials)
                norm = np.sqrt(pos.dot(pos.T))
                ground_truth[i,j,k] = np.sum(np.square(pos) / factor)
                gradient[i, j, k] = 2 * pos / factor
                gradient_abs[i, j, k] = np.linalg.norm(gradient[i, j, k])

    # gradient[np.isnan(gradient)] = 0.
    # gradient_abs[np.isnan(gradient_abs)] = 0.
    x = np.arange(tam) - initials[1]
    y = np.arange(tam) - initials[0]
    xx, yy = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.contour(ground_truth[0, :, :], xx, yy, 15, zdir='x', offset=-25)
    ax.contour(xx, ground_truth[:, 0, :], yy, 15, zdir='y', offset=25)
    ax.contour(xx, yy, ground_truth[:, :, 0], 15, zdir='z', offset=-25)
    ax.set_xlim3d(-25, 25)
    ax.set_ylim3d(-25, 25)
    ax.set_zlim3d(-25, 25)
    plt.title('$T_8$', fontsize=32)
    plt.show()

    inner_nodes_list = [1, 2]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
    for inner_nodes in inner_nodes_list:
        print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
        aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
