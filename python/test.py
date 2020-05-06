import mm
import numpy as np


if __name__ == '__main__':

    image = np.sqrt(3) / 3 * np.ones((50, 50, 50, 3))
    initials = [0, 0, 0]
    h = 1
    interp = 'pchip'
    search = 'gradient'

    u_surface = mm.wmm3D(image=image, initials=initials, h=h, interp=interp, search=search)

    print(u_surface[:5, :5, :5])



