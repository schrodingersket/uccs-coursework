import numpy as np


def compute_basic_solutions(A, b, level=0, start=0, basis=None, var_names=None):
    # Assert A and b are 2d and 1d arrays, respectively
    #
    assert(len(A.shape) == 2)
    assert(len(b.shape) == 1)

    A_rows, A_cols = A.shape

    assert(b.shape[0] == A_rows)

    if var_names:
        assert(len(var_names) == A_cols)

    if basis is None:
        basis = [None] * A_rows

    # Recursively compute all permutations of basic solutions
    #
    solutions = []
    for i in range(start, A_cols):
        basis[level] = i
        if level < A_rows - 1:
            solutions.extend(compute_basic_solutions(
                A,
                b,
                level=level + 1,
                start=(i + 1),
                basis=basis,
                var_names=var_names,
            ))
        else:
            B = np.array([A[:, beta] for beta in basis]).T
            # print('Basis: {}'.format([var_names[var] for var in basis] if var_names else basis))
            # print('{}'.format(B))
            try:
                x = np.linalg.solve(B, b)
                solutions.append({
                    'basis': basis.copy(),
                    'formatted_basis': [var_names[var] for var in basis] if var_names else basis,
                    'solution': x,
                })
                # print('Solution: {}'.format(x))

            except np.linalg.LinAlgError:
                solutions.append({
                    'basis': basis.copy(),
                    'formatted_basis': [var_names[var] for var in basis] if var_names else basis,
                    'solution': None,
                })
            #     print('Matrix is singular; two or more columns are linearly dependent.')
            # print('')

    return solutions


if __name__ == '__main__':
    # Test against Griva, Nash, Sofer Example 4.3
    #
    A = np.array((
        (-2, 1, 1, 0, 0),
        (-1, 1, 0, 1, 0),
        ( 1, 0, 0, 0, 1)
    ))

    b = np.array((
        2,
        3,
        3
    ))

    objective = lambda x1, x2, *_: -x1 - 2*x2
    minimum = None

    rows, cols = A.shape
    solns = compute_basic_solutions(A, b, var_names=['x1', 'x2', 's1', 's2', 's3'])

    for s in solns:
        soln = s.get('solution')
        pretty_basis = s.get('formatted_basis')
        print('Basis: {}'.format(pretty_basis))

        if soln is not None:
            full_solution = [0] * cols
            for i, beta in enumerate(s.get('basis')):
                full_solution[beta] = soln[i]

            feasible = all(i >= 0 for i in full_solution)
            print('Solution ({}feasible): {}'.format('' if feasible else 'in', full_solution))
            z = objective(*full_solution)
            print('Objective: {}'.format(z))

            if feasible and (minimum is None or z <= minimum[2]):
                if z == minimum[2]:
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
