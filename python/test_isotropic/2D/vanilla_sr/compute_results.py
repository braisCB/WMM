import numpy as np
import time
import mm
import cv2


def cart2pol(out):
    rho = np.sqrt(out[..., 0]**2 + out[..., 1]**2)
    phi = np.arctan2(out[..., 1], out[..., 0])
    return rho, np.cos(phi), np.sin(phi)

def pol2cart(rho, cos_phi, sin_phi):
    n = np.sqrt(cos_phi ** 2 + sin_phi ** 2)
    out = np.zeros(rho.shape + (2, ))
    out[..., 0] = rho * cos_phi / n
    out[..., 1] = rho * sin_phi / n
    return out

def aug_compute_results(gradient, ground_truth, initial_points, niter=10, inner_nodes=3, h=1):

    info = {}
    methods = ['fmm2D', 'msfm2D']
    orders = [1, 2]

    new_shape = tuple(((np.array(gradient.shape[:2]).astype(int) - 1) * inner_nodes + 1).tolist())
    new_initial_points = initial_points * inner_nodes
    new_h = h/inner_nodes

    rho, cos_phi, sin_phi = cart2pol(gradient)
    rho = cv2.resize(rho, dsize=new_shape, interpolation=cv2.INTER_CUBIC)
    cos_phi = cv2.resize(cos_phi, dsize=new_shape, interpolation=cv2.INTER_CUBIC)
    sin_phi = cv2.resize(sin_phi, dsize=new_shape, interpolation=cv2.INTER_CUBIC)
    new_gradient = pol2cart(rho, cos_phi, sin_phi)

    for method in methods:
        for order in orders:
            label = method + '_' + str(order)
            func = getattr(mm, method)
            print(label)
            info[label] = {}
            start_time = time.time()
            for _ in range(niter):
                info[label]['surface'] = func(new_gradient, new_initial_points, new_h, order)[::inner_nodes, ::inner_nodes]
            info[label]['time'] = (time.time() - start_time) / niter
            info[label]['metrics'] = {
                'l1': np.abs(info[label]['surface'] - ground_truth).mean(),
                'l2': np.square(info[label]['surface'] - ground_truth).mean(),
                'l_inf': np.abs(info[label]['surface'] - ground_truth).max(),
                'l_n_inf': (np.abs(info[label]['surface'] - ground_truth) / np.maximum(1e-10, ground_truth)).max()
            }

    for label, data in info.items():
        print(label, info[label]['metrics'], 'time : ', info[label]['time'])

