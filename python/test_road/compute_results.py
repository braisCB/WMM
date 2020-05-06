import numpy as np
import time
import distance


def fmm_compute_results(gradient, ground_truth, initial_points, niter, h=1,
                    fmm_techniques=('fmm2D', 'msfm2D'),
                    orders=(1, 2), factor=(1, 1)):

    info = {}
    factor = np.asarray(factor, dtype='int')

    final_points = np.array([k for k, _ in ground_truth], dtype=np.int32)
    initial_points = np.asarray(initial_points, dtype=np.int32) // factor
    final_points = np.asarray(final_points, dtype=np.int32) // factor
    gradient = gradient[::factor[0], ::factor[1], :]
    h = h * factor

    for fmm_technique in fmm_techniques:
        func = getattr(distance, fmm_technique)
        for order in orders:
            label = fmm_technique + '_' + str(order)
            print(label)
            info[label] = {}
            start_time = time.time()
            for _ in range(niter):
                info[label]['time'], info[label]['distance'] = func(gradient, initial_points, final_points, h, order)
            info[label]['cpu_time'] = (time.time() - start_time) / niter
            error = [info[label]['distance'][k[0] // factor[0], k[1] // factor[1]] - np.asarray(v) for k,v in ground_truth]
            info[label]['metrics'] = {
                'error': error,
                'l1': np.abs(error).mean(),
                'l2': np.square(error).mean(),
                'l_inf': np.abs(error).max()
            }
            print(label, info[label]['metrics'], 'time : ', info[label]['cpu_time'])

    return info


def compute_results(gradient, ground_truth, initial_points, niter, h=1,
                    interpolation_techniques=('linear', 'quad', 'pchip', 'hermite', 'spline'),
                    search_modes=('hopf_lax', 'golden_search', ), factor=1):

    info = {}
    factor = np.asarray(factor, dtype='int')

    final_points = np.array([k for k, _ in ground_truth], dtype=np.int32)
    initial_points = np.asarray(initial_points, dtype=np.int32) // factor
    final_points = np.asarray(final_points, dtype=np.int32) // factor
    gradient = gradient[::factor[0], ::factor[1], :]
    h *= factor

    for interpolation_technique in interpolation_techniques:
        for search_mode in search_modes:
            label = 'wmm_' + interpolation_technique + '_' + search_mode
            print(label)
            info[label] = {}
            start_time = time.time()
            for _ in range(niter):
                info[label]['time'], info[label]['distance'] = distance.wmm2D(gradient, initial_points, final_points, h, interpolation_technique, search_mode)
            info[label]['cpu_time'] = (time.time() - start_time) / niter
            error = [info[label]['distance'][k[0] // factor[0], k[1] // factor[1]] - np.asarray(v) for k,v in ground_truth]
            info[label]['metrics'] = {
                'error': error,
                'l1': np.abs(error).mean(),
                'l2': np.square(error).mean(),
                'l_inf': np.abs(error).max()
            }
            print(label, info[label]['metrics'], 'time : ', info[label]['cpu_time'])

    return info
