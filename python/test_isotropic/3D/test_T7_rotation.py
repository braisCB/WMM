import numpy as np
from compute_results import compute_results, aug_compute_results


if __name__ == '__main__':

    h = 1.
    tam = 51
    degrees = [0, 21, 39, 57, 75]

    initials = tam//2 * np.ones(3, dtype=int)

    for gamma in degrees:
        for beta in degrees:

            gamma_radian = gamma * np.pi / 180.
            beta_radian = beta * np.pi / 180.
            gamma_rotation = np.array([[1., 0., 0.], [0., np.cos(gamma_radian), -np.sin(gamma_radian)], [0., np.sin(gamma_radian), np.cos(gamma_radian)]])
            beta_rotation = np.array([[np.cos(beta_radian), -np.sin(beta_radian), 0.], [np.sin(beta_radian), np.cos(beta_radian), 0.], [0., 0., 1.]])

            gradient = np.empty((tam, tam, tam, 3))
            ground_truth = np.empty((tam, tam, tam))
            gradient_abs = np.empty((tam, tam, tam))

            factor = np.array([25, 16, 36])

            for i in range(tam):
                for j in range(tam):
                    for k in range(tam):
                        pos = h * (np.array([i, j, k]) - initials)
                        pos = beta_rotation @ gamma_rotation @ pos
                        norm = np.sqrt(pos.dot(pos.T))
                        ground_truth[i,j,k] = (np.square(pos) / factor).sum()
                        gradient[i, j, k] = gamma_rotation.T @ beta_rotation.T @ (2 * pos / factor)
                        gradient_abs[i, j, k] = np.linalg.norm(gradient[i, j, k])

            # gradient[np.isnan(gradient)] = 0.
            # gradient_abs[np.isnan(gradient_abs)] = 0.

            inner_nodes_list = [1, 2]

            niter = 1

            print('GAMMA ROTATION : ', gamma, 'DEGREES, BETA ROTATION : ', beta, ' DEGREES')
            compute_results(gradient=gradient, ground_truth=ground_truth,
                            initial_points=initials, niter=niter, h=h,
                            search_modes = ('gradient', ))
            # for inner_nodes in inner_nodes_list:
            #     print('Computing SR-WMM. Inner Nodes = ', inner_nodes)
            #     aug_compute_results(gradient=gradient, ground_truth=ground_truth,
            #                         initial_points=initials, niter=niter, inner_nodes=inner_nodes, h=h)
            # aug_compute_results(gradient=gradient, gradient_abs=gradient_abs, ground_truth=ground_truth,
            #                 initial_points=initials, niter=niter, nvalues=radius, h=h)
