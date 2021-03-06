import sys
sys.path.append('../../../')
import numpy as np
from compute_results import compute_results
from matplotlib import pyplot as plt


if __name__ == '__main__':

    h = 1
    tam = 101

    initials = tam//2 * np.ones((1,2), dtype=int)

    gradient = np.empty((tam, tam, 2))
    ground_truth = np.empty((tam, tam))
    gradient_abs = np.empty((tam, tam))

    x = np.arange(tam) - initials[0, 1]
    y = np.arange(tam) - initials[0, 0]
    xx, yy = np.meshgrid(x, y)

    for i in range(tam):
        for j in range(tam):
            diff = h * (np.array([i, j]) - initials)
            norm = np.sqrt(diff.dot(diff.T))
            ground_truth[i,j] = norm
            if norm == 0:
                diff = np.array([1., 0.])
                norm = 1
            gradient[i, j] = diff / norm
            gradient_abs[i, j] = np.linalg.norm(gradient[i, j])

    plt.figure()
    plt.contour(xx, yy, ground_truth, 30)
    plt.title('$T_1$', fontsize=32)
    plt.show()

    inner_nodes_list = [3, 5]

    niter = 1

    print('Computing WMM')
    compute_results(gradient=gradient, ground_truth=ground_truth,
                    initial_points=initials, niter=niter, h=h)
