import sys

sys.path.append('../')
import numpy as np
import imageio
from matplotlib import pyplot as plt
import compute_results as cr


if __name__ == '__main__':
    filename = './maps/galicia.png'
    filename_120 = './maps/galicia_120.png'
    filename_80 = './maps/galicia_80.png'
    filename_50 = './maps/galicia_50.png'
    output = '~/Documents/highway.png'

    image = np.asarray(imageio.imread(filename))[...,:3]
    mask_120 = np.expand_dims(np.asarray(imageio.imread(filename_120)), axis=2)
    mask_80 = np.expand_dims(np.asarray(imageio.imread(filename_80)), axis=2)
    mask_50 = np.expand_dims(np.asarray(imageio.imread(filename_50)), axis=2)

    # r, g, b = np.split(image, 3, axis=2)
    # r[mask_50 > 0] = 0
    # g[mask_50 > 0] = 0
    # b[mask_50 > 0] = 0
    # r[mask_80 > 0] = 0
    # g[mask_80 > 0] = 0
    # b[mask_80 > 0] = 255
    # r[mask_120 > 0] = 255
    # g[mask_120 > 0] = 0
    # b[mask_120 > 0] = 0
    # new_image = np.concatenate((r, g, b), axis=2).astype(np.uint8)
    # imageio.imwrite('./maps/galicia_roads.png', new_image, 'png')
    mask = 1000000. * np.ones_like(mask_120)
    mask[mask_50 > 0] = 60. / 40.
    mask[mask_80 > 200] = 60. / 80.
    mask[mask_120 > 200] = 60. / 120.

    corunha = (850, 1031)
    santiago = (1582, 892)
    vigo = (2488, 711)
    malpica = (884, 586)
    lugo = (1338, 1933)
    ourense = (2364, 1601)

    image_mask = np.zeros_like(mask, dtype=np.uint8)
    image_mask[mask_50 > 0] = 255
    image_mask[mask_80 > 200] = 255
    image_mask[mask_120 > 200] = 255

    h = 0.0794841368106578

    # mask = np.zeros(image.shape[:2] + (1,))
    # mask[(image[..., 0] > 250) & (image[..., 1] > 250) & (image[..., 2] > 250)] = 255
    # # #mask[image[..., 0] > 245] = 1.025 # 60./30.
    # # mask[(image[..., 0] > 245) & (image[...,2] < 200)] = 255 # 60./100.
    # mask = mask.astype(int)
    imageio.imwrite(output, image_mask[..., 0], 'png')

    mask = np.concatenate((mask, np.zeros_like(mask)), axis=2)

    # print(image.shape)
    # # ground_truth_gt = mm.ani_wmm2D(gradient_gt, initials_gt, h_gt, 'pchip', 'golden_search')
    # plt.figure()
    # plt.imshow(mask[...,0]/mask.max())
    # plt.show()

    ground_truth_distance = [(santiago, 71.1),
                             (malpica, 49.1),
                             (lugo, 91.),
                             (ourense, 170.),
                             (vigo, 153.)]

    ground_truth_time = [(santiago, 42.),
                         (lugo, 66.),
                         (ourense, 113.),
                         (vigo, 89.)]

    ground_truth = ground_truth_distance

    for factor in [(3, 3), (4, 4), (6, 6)]:
        print('FACTOR : ', factor)
        cr.fmm_compute_results(mask, ground_truth, [corunha], 1, h=h, factor=factor)
        cr.compute_results(mask, ground_truth, [corunha], 1, h=h, factor=factor)
