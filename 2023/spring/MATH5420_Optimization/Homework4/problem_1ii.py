import numpy as np

from simplex import simplex_algorithm


A = np.array((
  (4, 5,  2, 1, 1, 0),
  (3, 4, -1, 1, 0, 1),
))

b = np.array((
  20,
  30,
))

c = np.array((
  -5,
  -3,
  -2,
   0,
   0,
   0,
))

init = np.array((
  4, 
  5,
))

optimal_basis = simplex_algorithm(A, b, c, init)

if optimal_basis:
    xNi, xBi, xB, z = optimal_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\nOptimal basic feasible solution:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('optimal point = {}'.format(optimal_point))
    print('z = {}'.format(z))