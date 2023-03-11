import numpy as np

from simplex import simplex_algorithm

# Phase I
#
A = np.array((
  (2, -1, 3, 1, 0),
  (1,  2, 4, 0, 1),
))

b = np.array((
  30,
  40,
))

c = np.array((
  0,
  0,
  0,
  0,
  1,
))

init = np.array((
  4,
  3,
))

feasible_basis = simplex_algorithm(A, b, c, init, artificial_variables=[4])

if feasible_basis:
    xNi, xBi, xB, z = feasible_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\nBasic feasible solution for phase II:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('feasible point = {}'.format(optimal_point))
    print('z = {}'.format(z))

    print('\n\nEntering simplex phase II with xb = {}:'.format(xBi + 1))
    A = np.array((
        (2, -1, 3, 1),
        (1,  2, 4, 0),
    ))

    b = np.array((
        30,
        40,
    ))

    # Phase I
    #
    c = np.array((
        -4,
        -2,
        -8,
         0,
    ))

    optimal_basis = simplex_algorithm(A, b, c, xBi, artificial_variables=[])
    xNi, xBi, xB, z = optimal_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\nOptimal solution for phase II:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('optimal point = {}'.format(optimal_point))
    print('z = {}'.format(z))
