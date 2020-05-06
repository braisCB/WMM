import sys
sys.path.append('../../')
import numpy as np
import mm
from matplotlib import pyplot as plt


# def get_path(init, end, contour):
#     dirs = np.array([(-1, 0), (1, 0), (0, -1), (0, 1)])
#     diag_dirs = np.array([(-1, -1), (-1, 1), (1, -1), (1, 1)])
#     c_dirs = np.array([(0, 0), (1, 0), (0, 1), (1, 1)])
#
#     def contains(p):
#         contour_shape = np.array(contour.shape) - 1
#         return np.all(p >= 0, axis=-1) * np.all(p <= contour_shape, axis=-1)
#
#     def get_c(p):
#         p_floor = np.floor(p).astype(int)
#         diff = p - p_floor
#         p_neighs = p_floor + c_dirs
#         valid_p = contains(p_neighs)
#         c_neighs = np.zeros((4,))
#         c_neighs[valid_p] = contour[tuple(p_neighs[valid_p].T)]
#         if not valid_p[1]:
#             diff[0] = 0
#         if not valid_p[2]:
#             diff[1] = 0
#         c_p = (1. - diff[0])*(1. - diff[1])*c_neighs[0] + \
#               diff[0]*(1. - diff[1])*c_neighs[1] + \
#               (1. - diff[0])*diff[1]*c_neighs[2] + \
#               diff[0]*diff[1]*c_neighs[3]
#         return c_p
#
#     def get_gradient(p, c_p):
#         p_neighs = np.round(dirs + p).astype(int)
#         diff = p_neighs - p
#         wider = np.any(np.abs(diff) > 1, axis=-1)
#         if np.any(wider):
#             p_neighs[wider] -= dirs[wider]
#             diff[wider] = p_neighs[wider] - p
#         valid_p = contains(p_neighs)
#         c_neighs = np.Inf * np.ones((4,))
#         c_neighs[valid_p] = contour[tuple(p_neighs[valid_p].T)]
#         pos_y, pos_x = np.argmin(c_neighs[:2]), np.argmin(c_neighs[2:])
#         g_y = max(0., c_p - c_neighs[pos_y])
#         g_x = max(0., c_p - c_neighs[2 + pos_x])
#         if g_y == 0. and g_x == 0.:
#             p_neighs = np.round(diag_dirs + p).astype(int)
#             diff = p_neighs - p
#             wider = np.any(np.abs(diff) > 1, axis=-1)
#             if np.any(wider):
#                 p_neighs[wider] -= dirs[wider]
#                 diff[wider] = p_neighs[wider] - p
#             valid_p = contains(p_neighs)
#             c_neighs = np.Inf * np.ones((4,))
#             c_neighs[valid_p] = contour[tuple(p_neighs[valid_p].T)]
#             pos = np.argmin(c_neighs)
#             return diff[pos]
#         return np.array([g_y / diff[pos_y][0], g_x / diff[2 + pos_x][1]])
#
#     def normalize(p, grad):
#         grad_abs = np.abs(grad)
#         diff = p - np.round(p)
#         i = 0 if grad_abs[0] > grad_abs[1] else 1
#         factor = grad_abs[i] / np.abs(diff[i]) if grad[i] * diff[i] < 0 else grad_abs[i] / (1. - np.abs(diff[i]))
#         return factor
#
#     p = end.astype(float)
#     p_int = end
#     c_p = contour[tuple(p_int)]
#     path = [p_int]
#
#     while np.sum(np.abs(p_int - init) > 1):
#         print(p_int)
#         grad = get_gradient(p, c_p)
#         grad_1 = grad / normalize(p, grad)
#         p_aux = p + grad_1
#         c_p_aux = get_c(p_aux)
#         grad_aux = get_gradient(p_aux, c_p_aux)
#         grad_2 = 0.5 * (grad_aux + grad)
#         # grad_2[grad * grad_aux < 0] = 0
#         grad_2 /= normalize(p, grad_2)
#         p += grad_2
#         p_int = np.round(p).astype(int)
#         c_p = get_c(p)
#         if np.any(p_int != path[-1]):
#             path.append(p_int)
#     return np.array(path)


def get_path(init, end, contour, h):
    dirs = np.array([(-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1)], dtype=int)
    dirs_norm = np.linalg.norm(dirs, axis=-1)

    def contains(p):
        contour_shape = np.array(contour.shape) - 1
        return np.all(p >= 0, axis=-1) * np.all(p <= contour_shape, axis=-1)

    p = end
    path = [h * (p - initials)]
    c_neighs = np.zeros((8,))
    c_p = contour[tuple(p)]

    while np.any(p != init):
        p_neighs = p + dirs
        p_valid = contains(p_neighs)
        c_neighs[~p_valid] = np.Inf
        c_neighs[p_valid] = contour[tuple(p_neighs[p_valid].T)]
        grad = (c_neighs - c_p) / dirs_norm
        min_pos = np.argmin(grad)
        p, c_p = p_neighs[min_pos], c_neighs[min_pos]
        if np.any(p != path[-1]):
            path.append(h * (p - initials))
    return np.array(path)


if __name__ == '__main__':

    tam = 601
    h = 2 / (tam - 1)

    initials = tam//2 * np.ones(2, dtype=int)

    gradient = np.zeros((tam, tam, 4))
    gradient[:, :, 0] = 1
    gradient[:, :, 3] = 1

    w_0 = 6. * np.pi
    r_0 = .01
    i_d_0 = 5.
    d_0 = .01

    t_n = 4000

    t = np.linspace(0, 1, t_n)
    tubular = np.array([t * np.sin(w_0 * t), t * np.cos(w_0 * t)]).T

    for y in range(tam):
        for x in range(tam):
            pos = h * (np.array([y, x]) - initials)
            diff = np.linalg.norm(tubular - pos, axis=-1)
            mpos = np.argmin(diff)
            if diff[mpos] < r_0:
                g_y = np.sin(w_0 * t[mpos]) + t[mpos] * w_0 * np.cos(w_0 * t[mpos])
                g_x = np.cos(w_0 * t[mpos]) - t[mpos] * w_0 * np.sin(w_0 * t[mpos])
                g = np.array([g_y, g_x])
                g /= np.linalg.norm(g)
                D = np.array([[d_0*d_0, 0], [0, 1]])
                P = np.array([[g[0], -g[1]], [g[1], g[0]]])
                A = P @ D @ P.T
                gradient[y, x, 0] = A[0, 0]
                gradient[y, x, 1] = A[0, 1]
                gradient[y, x, 2] = A[1, 0]
                gradient[y, x, 3] = A[1, 1]
                # gradient[y, x, 0] = (i_d_0 * (np.sin(w_0 * t[mpos]) + t[mpos] * w_0 * np.cos(w_0 * t[mpos])))
                # gradient[y, x, 1] = (i_d_0 * (np.cos(w_0 * t[mpos]) - t[mpos] * w_0 * np.sin(w_0 * t[mpos])))
                # gradient[y, x, 0] = gradient[y, x, 0] if gradient[y, x, 0] == 0 else 1 / gradient[y, x, 0]
                # gradient[y, x, 1] = gradient[y, x, 1] if gradient[y, x, 1] == 0 else 1 / gradient[y, x, 1]
    x = np.linspace(-1, 1, tam)
    y = np.linspace(-1, 1, tam)

    xx, yy = np.meshgrid(x, y)

    surface = mm.ani_wmm2D_riemann(gradient, initials, h, 'pchip', 'golden_search')
    # aug_surface_3 = mm.aug_ani_wmm2D_riemann(gradient, initials, h, 3, 'pchip', 'golden_search')
    # aug_surface_10 = mm.aug_ani_wmm2D_riemann(gradient, initials, h, 10, 'pchip', 'golden_search')
    end = np.array([initials[0], tam - 1])
    path = get_path(initials, end, surface, h)
    plt.figure()
    plt.plot(path[:, 1], path[:, 0], 'r')
    plt.contourf(xx, yy, surface, 30)
    plt.title('$601 \\times 601$', fontsize=32)
    # plt.figure()
    # plt.contourf(xx, yy, aug_surface_3, 30)
    # plt.title('SR-WMM$^3$', fontsize=32)
    # plt.figure()
    # plt.contourf(xx, yy, aug_surface_10, 30)
    # plt.title('SR-WMM$^{10}$', fontsize=32)
    plt.show()


