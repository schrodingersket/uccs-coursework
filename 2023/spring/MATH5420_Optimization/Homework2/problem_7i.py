import numpy as np

from vector_projection import proj

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

A = np.array((
    (1,  1,  1, 1),
    (1, -1, -1, 1),
    (0,  1,  0, 1),
))

x1 = np.array((1,  3,  1, 2))
x2 = np.array((0, -2, -3, 4))

bN = np.array((
    (-1, -1, 1, 1),
))

bR = np.array((
    (1,  1,  1, 1),
    (1, -1, -1, 1),
    (0,  1,  0, 1),
))

print('x1: {}'.format(x1))
c1, p1 = proj(x1, bN[0,:])
q1 = x1 - (c1*p1)

print('p1 = {} * {} = {}'.format(c1, p1, c1*p1))
print('q1 = x1 - p1 = {}'.format(q1))
print('')
print('x2: {}'.format(x2))
c2, p2 = proj(x2, bN[0,:])
q2 = x2 - (c2*p2)

print('p2 = {} * {} = {}'.format(c2, p2, c2*p2))
print('q2 = x2 - p2 = {}'.format(q2))