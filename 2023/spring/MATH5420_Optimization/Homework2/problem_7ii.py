import numpy as np

from vector_projection import proj

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

A = np.array((
    (1,  1,  1, 1),
))

x1 = np.array((-2, 4, 5, -2))
x2 = np.array((7, 5, -13, 1))

bN = np.array((
    (-1, 1, 0, 0),
    (-1, 0, 1, 0),
    (-1, 0, 0, 1),
))

bR = np.array((
    (1,  1,  1, 1),
))

print('x1: {}'.format(x1))
c1, q1 = proj(x1, bR[0,:])
p1 = x1 - (c1*q1)

print('q1 = {} * {} = {}'.format(c1, q1, c1*q1))
print('p1 = x1 - q1 = {}'.format(p1))
print('')
print('x2: {}'.format(x2))
c2, q2 = proj(x2, bR[0,:])
p2 = x2 - (c2*q2)

print('q2 = {} * {} = {}'.format(c2, q2, c2*q2))
print('p2 = x2 - q2 = {}'.format(p2))