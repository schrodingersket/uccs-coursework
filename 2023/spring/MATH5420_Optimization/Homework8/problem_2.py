import numpy as np

Z1 = np.array((
    1,
    2
))

DDF1 = np.array((
    (-2, 0),
    (0, 1),
))

print('a = 2 (reduced Hessian): {}'.format(np.matmul(Z1.T, np.matmul(DDF1, Z1))))

Z2 = np.array((
    2,
    1
))

DDF2 = np.array((
    (-1/2, 0),
    (0,    1),
))

print('a = -1 (reduced Hessian): {}'.format(np.matmul(Z2.T, np.matmul(DDF2, Z2))))
