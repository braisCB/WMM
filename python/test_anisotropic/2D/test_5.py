import sys
sys.path.append('../../')
import numpy as np
from compute_results_riemann import compute_results, aug_compute_results
import mm
from matplotlib import pyplot as plt
import pickle
import os


if __name__ == '__main__':

    tam_gt = 292 * 14 + 1
    inner_nodes_list = [3, 5]
    filename = 'test_5_results.pickle'
    gradient_gt = np.zeros((tam_gt, tam_gt, 4))
    gradient_gt_aux = np.zeros((tam_gt, tam_gt, 2))
    gradient_gt_aux_2 = np.zeros((tam_gt, tam_gt, 2))
    ps = np.zeros((tam_gt, tam_gt, 2))
    initials_gt = tam_gt // 2 * np.ones((1, 2), dtype=int)
    h_gt = 1 / (tam_gt - 1)

    degrees = np.arange(46, 46, dtype=int)
    xs = np.linspace(-0.5, 0.5, tam_gt)
    ps_aux = np.array(np.meshgrid(xs, xs)).transpose([1, 2, 0])

    methods = ['wmm']
    for inner_node in inner_nodes_list:
        methods.append('sr-wmm-' + str(inner_node))

    if os.path.exists(filename):
        print('loading file')
        with open(filename, 'rb') as handle:
            stats = pickle.load(handle)
    else:
        stats = {
            'angle': np.zeros_like(degrees),
        }
        for method in methods:
            stats[method] = {
                'l1': {
                    'gradient': np.Inf * np.ones_like(degrees),
                    'hopf_lax': np.Inf * np.ones_like(degrees),
                    'golden_search': np.Inf * np.ones_like(degrees),
                },
                'l_inf': {
                    'gradient': np.Inf * np.ones_like(degrees),
                    'hopf_lax': np.Inf * np.ones_like(degrees),
                    'golden_search': np.Inf * np.ones_like(degrees),
                }
            }

    for i, degree in enumerate(degrees):
        i = degree

        radian = degree * np.pi / 180.
        stats['angle'][i] = radian
        rotation_matrix = np.array([[np.cos(radian), np.sin(radian)], [-np.sin(radian), np.cos(radian)]])

        ps[:, :, 0] = rotation_matrix[0, 0] * ps_aux[:, :, 1] + rotation_matrix[0, 1] * ps_aux[:, :, 0]
        ps[:, :, 1] = rotation_matrix[1, 0] * ps_aux[:, :, 1] + rotation_matrix[1, 1] * ps_aux[:, :, 0]

        gradient_gt_aux[:, :, 0] = 0.75 * 3 * np.pi * np.cos(3 * np.pi * ps[:, :, 0]) * np.sin(3 * np.pi * ps[:, :, 1])
        gradient_gt_aux[:, :, 1] = 0.75 * 3 * np.pi * np.sin(3 * np.pi * ps[:, :, 0]) * np.cos(3 * np.pi * ps[:, :, 1])
        gradient_gt_aux_2[:, :, 0] = rotation_matrix[0, 0] * gradient_gt_aux[:, :, 0] + rotation_matrix[1, 0] * gradient_gt_aux[:, :, 1]
        gradient_gt_aux_2[:, :, 1] = rotation_matrix[0, 1] * gradient_gt_aux[:, :, 0] + rotation_matrix[1, 1] * gradient_gt_aux[:, :, 1]
        gradient_gt[:, :, 0] = 1. + gradient_gt_aux_2[:, :, 0] * gradient_gt_aux_2[:, :, 0]
        gradient_gt[:, :, 1] = gradient_gt[:, :, 2] = gradient_gt_aux_2[:, :, 0] * gradient_gt_aux_2[:, :, 1]
        gradient_gt[:, :, 3] = 1. + gradient_gt_aux_2[:, :, 1] * gradient_gt_aux_2[:, :, 1]
        # matrix = np.array([[gradient_gt[:, :, 0], gradient_gt[:, :, 1]], [gradient_gt[:, :, 2], gradient_gt[:, :, 3]]])
        # matrix = np.transpose(matrix, (2,3,0,1))
        # w, v = np.linalg.eigh(matrix)
        # gradient_gt[:, :, 0] = v[:, :, 1, 0]
        # gradient_gt[:, :, 1] = v[:, :, 1, 1]
        # gradient_gt[:, :, 2] = w[:, :, 1]
        # gradient_gt[:, :, 3] = w[:, :, 0]

        # for y in range(tam_gt):
        #     for x in range(tam_gt):
        #         p = (h_gt * (np.array([y, x]) - initials_gt))[0]
        #         p = rotation_matrix @ p
        #         gradient_gt2[y, x, 0] = 0.75 * 3 * np.pi * np.cos(3 * np.pi * p[0]) * np.sin(3 * np.pi * p[1])
        #         gradient_gt2[y, x, 1] = 0.75 * 3 * np.pi * np.sin(3 * np.pi * p[0]) * np.cos(3 * np.pi * p[1])
        #         gradient_gt2[y, x] = rotation_matrix.T @ gradient_gt2[y, x]
        print('Computing ground truth')
        ground_truth_gt = mm.ani_wmm2D_riemann(gradient_gt, initials_gt, h_gt, 'linear', 'golden_search')
        # ground_truth_gt = mm.oum2D(gradient_gt, initials_gt, h_gt, 5, 'golden_search')
        plt.contour(ground_truth_gt.T, 30)
        # plt.axis('off')
        plt.show()

        reductions = [14]
        for reduction in reductions:
            tam = (tam_gt - 1) / reduction + 1
            initials = tam // 2 * np.ones((1, 2), dtype=int)
            h = reduction * h_gt

            gradient = gradient_gt[::reduction, ::reduction, :]

            ground_truth = ground_truth_gt[::reduction, ::reduction]

            niter = 1
            print('SIZE : ', tam)
            print('DEGREE : ', degree)
            print('Computing WMM')
            info = compute_results(gradient=gradient, ground_truth=ground_truth,
                                   initial_points=initials, niter=niter, h=h)
            for label, data in info.items():
                key = '_'.join(label.split('_')[2:])
                stats['wmm']['l1'][key][i] = min(stats['wmm']['l1'][key][i], info[label]['metrics']['l1'])
                stats['wmm']['l_inf'][key][i] = min(stats['wmm']['l_inf'][key][i], info[label]['metrics']['l_inf'])
            for inner_node in inner_nodes_list:
                method = 'sr-wmm-' + str(inner_node)
                print('Computing SR-WMM. Inner Nodes = ', inner_node)
                info = aug_compute_results(gradient=gradient, ground_truth=ground_truth,
                                           initial_points=initials, niter=niter, inner_nodes=inner_node, h=h)
                for label, data in info.items():
                    key = '_'.join(label.split('_')[3:])
                    stats[method]['l1'][key][i] = min(stats[method]['l1'][key][i], info[label]['metrics']['l1'])
                    stats[method]['l_inf'][key][i] = min(stats[method]['l_inf'][key][i], info[label]['metrics']['l_inf'])

        with open(filename, 'wb') as handle:
            pickle.dump(stats, handle, protocol=pickle.HIGHEST_PROTOCOL)
        del ground_truth_gt

    legend = []
    stats['angle'] = np.arange(46) * np.pi / 180.
    for error_type in ['l1', 'l_inf']:
        error_label = '$L_1$' if error_type == 'l1' else '$L_{\\infty}$'
        plt.figure()
        # plt.plot(stats['angle'], stats['wmm'][error_type]['gradient'])
        plt.plot(stats['angle'], stats['wmm'][error_type]['hopf_lax'])
        plt.plot(stats['angle'], stats['wmm'][error_type]['golden_search'])
        legend += ['wmm$_{hl}$', 'wmm$_{gs}$']
        # legend += ['wmm$_{gs}$']
        for inner_node in [3, 5]:
            method = 'sr-wmm-' + str(inner_node)
            key = 'sr-wmm$^' + str(inner_node)
            # plt.plot(stats['angle'], stats[method][error_type]['gradient'])
            plt.plot(stats['angle'], stats[method][error_type]['hopf_lax'])
            plt.plot(stats['angle'], stats[method][error_type]['golden_search'])
            legend += [key + '_{hl}$', key + '_{gs}$']
            # legend += [key + '_{gs}$']
        plt.title(error_label + ' error', fontsize=32)
        plt.xlabel('$\\theta$', fontsize=22)
        plt.legend(legend, fontsize=12)
    plt.show()
