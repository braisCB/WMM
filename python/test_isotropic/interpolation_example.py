from scipy.interpolate import PchipInterpolator, interp1d
from matplotlib import pyplot as plt
import numpy as np


if __name__ == '__main__':
    x = np.array([0, 1, 2])
    y = np.array([0, 1, 10])

    keys = ['linear', 'cubic', 'pchip']

    plt.figure()
    plt.plot(x, y, '*')
    cubic = interp1d(x, y, kind='cubic')
    x_new = np.linspace(0., 2., 1000)
    y_cubic = cubic(x_new)
    pchip = PchipInterpolator(x, y)
    y_pchip = pchip(x_new)
    plt.plot(x_new, y_cubic)
    plt.plot(x_new, y_pchip)
    plt.title('Interpolation Problem')
    plt.legend(keys)
    plt.show()

