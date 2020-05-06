import sys
sys.path.append('../../../')
import numpy as np
from compute_results import compute_results
from matplotlib import pyplot as plt


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

    x = np.arange(tam) - initials[1]
    y = np.arange(tam) - initials[0]
    xx, yy = np.meshgrid(x, y)
    plt.figure()
    plt.contour(xx, yy, ground_truth, 30)
    plt.title('$T_3$', fontsize=32)
    plt.show()
    # gradient[np.isnan(gradient)] = 0.
    # gradient_abs[np.isnan(gradient_abs)] = 0.

    inner_nodes_list = [3, 5]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
