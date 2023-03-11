import numpy as np

from simplex import simplex_algorithm

A = np.array((
  ( 2, -3, 2, -1,  0, 1, 0),
  (-1,  1, 1,  0, -1, 0, 1),
))

b = np.array((
  3,
  5,
))

# Phase I
#
c = np.array((
  0,
  0,
  0,
  0,
  0,
  1,
  1,
))

init = np.array((
  5,
  6,
))

feasible_basis = simplex_algorithm(A, b, c, init, artificial_variables=[4])

if feasible_basis:
    xNi, xBi, xB, z = feasible_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\nBasic feasible solution for phase I:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('feasible point = {}'.format(optimal_point))
    print('z = {}'.format(z))
