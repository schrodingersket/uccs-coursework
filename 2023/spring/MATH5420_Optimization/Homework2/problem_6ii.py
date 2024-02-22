import numpy as np

from ratio_test import ratio_test

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

A = [
    np.array((9, 4, 1, 9, -7)),
    np.array((6, -7, 8, -4, -6)),
    np.array((1, 6, 3, -7, 6)),
]

b = [
    np.array((-15,)),
    np.array((-30,)),
    np.array((-20,)),
]

x = np.array((7, -4, -3, -3, 3))
p = np.array((3, 2, 0, 1, -2))


ratio_test(A, b, x, p)