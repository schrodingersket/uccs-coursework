import numpy as np

from basic_solutions import compute_basic_solutions

A = np.array((
    (-3, 2, 1, 0),
    (-2, 1, 0, 1)
))

b = np.array((
    30,
    12,
))

objective = lambda x1, x2, *_: -5*x1 - 7*x2
minimum = None

rows, cols = A.shape
solns = compute_basic_solutions(A, b, var_names=['x1', 'x2', 's1', 's2'])

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
        z = objective(*full_solution)
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
