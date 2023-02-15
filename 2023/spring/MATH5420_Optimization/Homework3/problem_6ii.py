import numpy as np

from basic_solutions import compute_basic_solutions

# Test against Griva, Nash, Sofer Example 4.3
#
A = np.array((
    (2, 1, 1, 0, 0, 0),
    (3, 1, 0, 1, 0, 0),
    (4, 1, 0, 0, 1, 0),
    (5, 1, 0, 0, 0, 1)
))

b = np.array((
    3,
    4,
    5,
    6
))

rows, cols = A.shape
solns = compute_basic_solutions(A, b, var_names=['x1', 'x2', 's1', 's2', 's3', 's4'])

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
    else:
        print('No solution exists!')

    print('')
    input('Press [Enter] to see the next point:')