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
    for i in range(start, A_cols):
        basis[level] = i
        if level < A_rows - 1:
            compute_basic_solutions(
                A,
                b,
                level=level + 1,
                start=(i + 1),
                basis=basis,
                var_names=var_names,
            )
        else:
            B = np.array([A[:, beta] for beta in basis]).T
            print('Basis: {}'.format([var_names[var] for var in basis] if var_names else basis))
            print('{}'.format(B))
            try:
                print('Solution: {}'.format(np.linalg.solve(B, b)))
            except np.linalg.LinAlgError:
                print('Matrix is singular; two or more columns are linearly dependent.')
            print('')


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

    compute_basic_solutions(A, b, var_names=['x1', 'x2', 's1', 's2', 's3'])
