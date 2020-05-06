import sys
sys.path.append('../../')
import numpy as np
import mm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_path(init, end, propagation, h, zero_point):
    p = end
    path = [h * (p - zero_point)]
    while np.any(p != init):
        dir = propagation[tuple(p)]
        if np.abs(dir).sum() == 0:
            break
        p -= dir
        path.append(h * (p - zero_point))
    return np.array(path)

# def get_path(init, end, contour, h, zero_point):
#     dirs = []
#     for i in range(-1, 2):
#         for j in range(-1, 2):
#             for k in range(-1, 2):
#                 if i == 0 and j == 0 and k == 0:
#                     continue
#                 dirs.append((i, j, k))
#     dirs = np.array(dirs)
#     dirs_norm = np.linalg.norm(dirs, axis=-1)
#
#     def contains(p):
#         contour_shape = np.array(contour.shape) - 1
#         return np.all(p >= 0, axis=-1) * np.all(p <= contour_shape, axis=-1)
#
#     p = end
#     path = [h * (p - zero_point)]
#     c_neighs = np.zeros((len(dirs),))
#     c_p = contour[tuple(p)]
#
#     while np.any(p != init):
#         p_neighs = p + dirs
#         p_valid = contains(p_neighs)
#         c_neighs[~p_valid] = np.Inf
#         c_neighs[p_valid] = contour[tuple(p_neighs[p_valid].T)]
#         grad = (c_neighs - c_p) / dirs_norm
#         min_pos = np.argmin(grad)
#         p, c_p = p_neighs[min_pos], c_neighs[min_pos]
#         if np.any(p != path[-1]):
#             path.append(h * (p - zero_point))
#     return np.array(path)


if __name__ == '__main__':

    tam = (201, 201, 273)
    h = 2.2 / (tam[0] - 1)

    initials = (191, 100, 0)
    zero_point = (100, 100, 0)

    gradient = np.zeros(tam + (9, ))
    gradient[:, :, :, 0] = 1
    gradient[:, :, :, 4] = 1
    gradient[:, :, :, 8] = 1

    w_0 = 5./2. * np.pi
    r_0 = .02
    i_d_0 = 5.
    d_0 = .02

    t_n = 8000

    t = np.linspace(-0.1, 3.1, t_n)
    tubular = np.array([np.cos(w_0 * t), np.sin(w_0 * t), t]).T
    cont = 0

    for z in range(tam[2]):
        print(z, cont)
        for x in range(tam[0]):
            for y in range(tam[1]):
                pos = h * (np.array([x, y, z]) - zero_point)
                diff = np.linalg.norm(tubular - pos, axis=-1)
                mpos = np.argmin(diff)
                if diff[mpos] < r_0:
                    g_x = - w_0 * np.sin(w_0 * t[mpos])
                    g_y = w_0 * np.cos(w_0 * t[mpos])
                    g_z = 1
                    g_0 = np.array([g_x, g_y, g_z])
                    g_1 = np.array([-g_y, g_x, 0])
                    g_2 = np.cross(g_0, g_1)
                    P = np.array((g_0, g_1, g_2)).T
                    P /= np.linalg.norm(P, axis=0, keepdims=True)
                    D = np.array([[d_0*d_0, 0, 0], [0, 1, 0], [0, 0, 1]])
                    # g = g_0 / np.linalg.norm(g_0)
                    # H = np.eye(3) - (2 * g @ g.T) / (g.T @ g)
                    # _, P = np.linalg.eigh(H)
                    # P[:, 0] = -1. * P[:, 0]
                    A = P @ D @ P.T
                    gradient[x, y, z, 0] = A[0, 0]
                    gradient[x, y, z, 1] = A[0, 1]
                    gradient[x, y, z, 2] = A[0, 2]
                    gradient[x, y, z, 3] = A[1, 0]
                    gradient[x, y, z, 4] = A[1, 1]
                    gradient[x, y, z, 5] = A[1, 2]
                    gradient[x, y, z, 6] = A[2, 0]
                    gradient[x, y, z, 7] = A[2, 1]
                    gradient[x, y, z, 8] = A[2, 2]
                    # gradient[x, y, z, 0] = (i_d_0 * (np.sin(w_0 * t[mpos]) + t[mpos] * w_0 * np.cos(w_0 * t[mpos])))
                    # gradient[x, y, z, 1] = (i_d_0 * (np.cos(w_0 * t[mpos]) - t[mpos] * w_0 * np.sin(w_0 * t[mpos])))
                    # gradient[x, y, z, 0] = gradient[x, y, z, 0] if gradient[x, y, z, 0] == 0 else 1 / gradient[x, y, z, 0]
                    # gradient[x, y, z, 1] = gradient[x, y, z, 1] if gradient[x, y, z, 1] == 0 else 1 / gradient[x, y, z, 1]
                    cont += 1

    print('Computing surface')
    surface, propagation = mm.ani_wmm3D_riemann(gradient, initials, h, 'pchip', 'golden_search')
    # aug_surface_3 = mm.aug_ani_wmm2D_riemann(gradient, initials, h, 3, 'pchip', 'golden_search')
    # aug_surface_10 = mm.aug_ani_wmm2D_riemann(gradient, initials, h, 10, 'pchip', 'golden_search')
    end = np.array([100, 9, 272])
    print('Computing path')
    path = get_path(initials, end, propagation, h, zero_point)
    print('Path computed')
    x = h * (np.arange(tam[0]) - initials[0])
    z = h * (np.arange(tam[2]) - initials[2])
    xx, zz = np.meshgrid(x, z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.contourf(surface[0, :, :], xx, zz, 30, zdir='x', offset=-1.1)
    # ax.contourf(xx, surface[:, 0, :], zz, 30, zdir='y', offset=1.1)
    # ax.contourf(xx, xx, surface[:, :, 0], 30, zdir='z', offset=0)

    ax.plot(path[:, 0], path[:, 1], path[:, 2], zdir='z')
    ax.set_xlim3d(-1.1, 1.1)
    ax.set_ylim3d(-1.1, 1.1)
    ax.set_zlim3d(0, 3)
    # plt.title('3D spiral', fontsize=32)
    plt.show()
