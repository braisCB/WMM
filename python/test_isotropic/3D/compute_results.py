import numpy as np
import time
import sys
sys.path.append('../../')
import mm


def compute_results(gradient, ground_truth, initial_points, niter=10, h=1,
                    interpolation_techniques = ('pchip', 'bilinear', 'spline', 'hermite', 'cubic'),
                    search_modes = ('gradient', 'hopf_lax', 'golden_search')):

    info = {}

    for interpolation_technique in interpolation_techniques:
        for search_mode in search_modes:
            label = 'wmm_' + interpolation_technique + '_' + search_mode
            print(label)
            info[label] = {}
            start_time = time.time()
            for _ in range(niter):
                info[label]['surface'] = mm.wmm3D(gradient, initial_points, h, interpolation_technique, search_mode)
            info[label]['time'] = (time.time() - start_time) / niter
            info[label]['metrics'] = {
                'l1': np.abs(info[label]['surface'] - ground_truth).mean(),
                'l2': np.square(info[label]['surface'] - ground_truth).mean(),
                'l_inf': np.abs(info[label]['surface'] - ground_truth).max(),
                'l_n_inf': (np.abs(info[label]['surface'] - ground_truth) / np.maximum(1e-10, ground_truth)).max()
            }
            print(label, info[label]['metrics'], 'time : ', info[label]['time'])

    for label, data in info.items():
        print(label, info[label]['metrics'], 'time : ', info[label]['time'])


def aug_compute_results(gradient, ground_truth, initial_points, niter=10, inner_nodes=3, h=1,
                        interpolation_techniques = ('pchip', 'bilinear', 'spline', 'hermite', 'cubic'),
                        search_modes = ('gradient', 'hopf_lax', 'golden_search')):

    info = {}

    for interpolation_technique in interpolation_techniques:
        for search_mode in search_modes:
            label = 'aug_wmm_' + interpolation_technique + '_' + search_mode
            print(label)
            info[label] = {}
            start_time = time.time()
            for _ in range(niter):
                info[label]['surface'] = mm.aug_wmm3D(gradient, initial_points, h, inner_nodes, interpolation_technique, search_mode)
            info[label]['time'] = (time.time() - start_time) / niter
            info[label]['metrics'] = {
                'l1': np.abs(info[label]['surface'] - ground_truth).mean(),
                'l2': np.square(info[label]['surface'] - ground_truth).mean(),
                'l_inf': np.abs(info[label]['surface'] - ground_truth).max(),
            }
            print(label, info[label]['metrics'], 'time : ', info[label]['time'])

    for label, data in info.items():
        print(label, info[label]['metrics'], 'time : ', info[label]['time'])

