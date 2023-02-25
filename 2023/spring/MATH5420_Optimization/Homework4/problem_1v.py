import numpy as np

from plot_feasible_region import plot_2d_region
from simplex import simplex_algorithm
from basic_solutions import compute_basic_solutions

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
    xlims=(-2, 100),
    ylims=(-2, 100),
    usetex=True,
    legend=True,
    savefile='problem_1v.png',
    objective_fn=lambda z, x1: (1/8)*z - (7/8)*x1,  # Solve z = 7*x_1 + 8*x_2 for x_2
)


optimal_basis = simplex_algorithm(A, b, c, init)

if optimal_basis:
    xB, xN, z = optimal_basis

    print('xb = {}'.format(xB + 1))
    print('xn = {}'.format(xN + 1))
    print('z = {}'.format(z))

objective = lambda x1, x2, *_: -7*x1 - 8*x2
minimum = None
