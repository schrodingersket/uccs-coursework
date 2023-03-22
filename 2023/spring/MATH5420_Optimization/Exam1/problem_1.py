import numpy as np

from basic_solutions import compute_basic_solutions
from simplex import simplex_algorithm

A = np.array((
  (1, 4, 2, 1, 0),
  (1, 2, 4, 0, 1),
))

b = np.array((
  48,
  60,
))

# Phase I
#
c = np.array((
  -6,
  -14,
  -13,
   0,
   0,
))

init = np.array((
  3,
  4,
))

feasible_basis = simplex_algorithm(A, b, c, init, verbose=True)

if feasible_basis:
    xNi, xBi, xB, z = feasible_basis
    optimal_point = np.zeros(c.shape)
    optimal_point[xBi] = xB

    print('\n\nBasic feasible solution via simplex:')
    print('xb = {}'.format(xBi + 1))
    print('xn = {}'.format(xNi + 1))
    print('feasible point = {}'.format(optimal_point))
    print('z = {}'.format(z))

print('\n\nBrute-forcing solutions...')


# Find basic solutions via brute force to verify simplex solution
#
rows, cols = A.shape
solns = compute_basic_solutions(A, b, var_names=['x1', 'x2', 'x3', 's4', 's5'])

minimum = None

for s in solns:
    soln = s.get('solution')
    pretty_basis = s.get('formatted_basis')
    print('Basis: {}'.format(pretty_basis))
    print('B: {}'.format(np.array([A[:, k] for k in s.get('basis')]).T))

    if soln is not None:
        full_solution = [0] * cols
        for i, beta in enumerate(s.get('basis')):
            full_solution[beta] = soln[i]

        feasible = all(i >= 0 for i in full_solution)
        print('Solution ({}feasible): {}'.format('' if feasible else 'in', full_solution))
        z = np.dot(full_solution, c)
        print('Objective: {}'.format(z))

        if feasible and (minimum is None or z <= minimum[2]):
            if minimum is not None and z == minimum[2]:
                print('[Degenerate solution]')
            else:
                minimum = (pretty_basis, full_solution, z)
    else:
        print('No solution exists!')

    print('')
    input('Press [Enter] to see the next point:')
    

print('Basic optimal feasible solution achieved at:')
beta, full_solution, z = minimum
print(beta)
print(full_solution)
print(z)