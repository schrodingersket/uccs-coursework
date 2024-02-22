import numpy as np

from plot_feasible_region import plot_2d_region
from simplex import simplex_algorithm

# Primal problem
#

# Phase I
#
print('\n\nSolving phase I primal problem: \n\n')
A = np.array((
    (-1, 1, -1, 0, 1),
    (2,  -1, 0, 1, 0),
))

b = np.array((
    1,
    2,
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
    print('\n\n[PRIMAL Phase I] Initial basic feasible solution:')
    print('z = {}'.format(z))
    print('Entering basis is xb = {}'.format(xBi + 1))
else:
    raise Exception('Problem is infeasible!')

print('\n\nSolving phase II primal problem: \n\n')
A = np.array((
  (1, -1, 1, 0),
  (2, -1, 0, 1),
))

b = np.array((
  -1,
   2,
))

c = np.array((
    1,
    1,
    0,
    0
))

init = xBi


plot_2d_region(
    A,
    b,
    c[0:2],
    xlims=(-2, 5),
    ylims=(-2, 5),
    usetex=True,
    legend=True,
    savefile='problem_5_primal.png',
)

optimal_basis = simplex_algorithm(A, b, c, init, verbose=True)

if optimal_basis:
    xNi, xBi, xB, z = optimal_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\n[PRIMAL Phase II] Optimal basic feasible solution:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('optimal point = {}'.format(optimal_point))
    print('-z = {}'.format(z))

# Dual problem
#
print('\n\nSolving dual problem: \n\n')
A = np.array((
    (-1, -2, 1, 0),
    ( 1,  1, 0, 1),
))

b = np.array((
    1,
    1,
))

c = np.array((
   -1,
    2,
    0,
    0,
    0
))

init = np.array((
    2,
    3,
))


plot_2d_region(
    A,
    b,
    c[0:2],
    xlims=(-2, 5),
    ylims=(-2, 5),
    usetex=True,
    legend=True,
    savefile='problem_5_dual.png',
    latex_vars=('y\'_1', 'y_2'),
    plaintext_vars=('y1\'', 'y2')
)

optimal_basis = simplex_algorithm(A, b, c, init)

if optimal_basis:
    xNi, xBi, xB, z = optimal_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\n[DUAL] Optimal basic feasible solution:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('optimal point = {}'.format(optimal_point))
    print('z = {}'.format(z))
