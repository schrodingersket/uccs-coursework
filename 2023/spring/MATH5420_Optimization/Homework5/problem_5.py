import numpy as np

from plot_feasible_region import plot_2d_region
from simplex import simplex_algorithm

A = np.array((
  (4, 1, 1, 0, 0),
  (1, 1, 0, 1, 0),
  (1, 0, 0, 0, 1),
))

b = np.array((
  100,
  80,
  40,
))

c = np.array((
    -7,
    -8,
     0,
     0,
     0,
))

init = np.array((
    2,
    3,
    4,
))


plot_2d_region(
    A,
    b,
    c[0:2],
    xlims=(-2, 100),
    ylims=(-2, 100),
    usetex=True,
    legend=True,
    savefile='problem_1v_1.png',
    active_point=(0, 0)
)

plot_2d_region(
    A,
    b,
    c[0:2],
    xlims=(-2, 100),
    ylims=(-2, 100),
    usetex=True,
    legend=True,
    savefile='problem_1v_2.png',
    active_point=(0, 80)
)


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
