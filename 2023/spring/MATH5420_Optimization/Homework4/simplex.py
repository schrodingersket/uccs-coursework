import numpy as np

class MaxIterationError(Exception):
    pass

def simplex_algorithm(A, b, c, x_b_idx, max_iters=10):
    optimal = False
    iters = 0
    while not optimal and iters < max_iters:
        x_n_idx = np.setdiff1d(np.arange(0, A.shape[1], 1), x_b_idx)

        B = A[:, x_b_idx]
        N = A[:, x_n_idx]
        c_b = c[x_b_idx]
        c_n = c[x_n_idx]

        x_b = np.linalg.solve(B, b)
        y = np.linalg.solve(B.T, c_b)

        c_n_hat = c_n.T - np.matmul(y.T, N)
        tableau_body = np.linalg.solve(B, N)
        z = np.matmul(np.matmul(c_b.T, np.linalg.inv(B)), b)

        entering_variable_idx = np.argmin(c_n_hat)
        entering_variable = x_n_idx[entering_variable_idx]
        entering_column = np.linalg.solve(B, A[:, entering_variable])

        print('B: {}'.format(B))
        print('N: {}'.format(N))
        print('c_b: {}'.format(c_b))
        print('c_n: {}'.format(c_n))
        print('c_n_hat (reduced cost): {}'.format(c_n_hat))
        # print('tableau: {}'.format(tableau_body))
        print('rhs (b_hat): {}'.format(x_b))
        # print('z: {}'.format(z))
        print('entering column (A_{}): {}'.format(entering_variable + 1, entering_column))
        # print('basis: {}'.format(x_b_idx + 1))

        if np.all(c_n_hat >= 0):
            optimal = True
            return x_n_idx, x_b_idx, x_b, z

        ratios = np.ma.array(x_b / entering_column, mask=(x_b / entering_column) <= 0)
        leaving_variable_idx = np.argmin(ratios)
        leaving_variable = x_b_idx[leaving_variable_idx]

        print('x{} is entering and x{} is leaving at iteration {}'.format(
            entering_variable + 1,
            leaving_variable + 1,
            iters + 1,
        ))

        x_b_idx[leaving_variable_idx] = entering_variable
        iters += 1
        if iters == max_iters:
            raise MaxIterationError('Exceeded max simplex iteration count! Is your problem bounded?')
        else:
            input('\n\nPress [Enter] for the next iteration: ')



if __name__ == '__main__':
    A = np.array((
        (-2, 1, 1, 0, 0),
        (-1, 2, 0, 1, 0),
        ( 1, 0, 0, 0, 1),
    ))

    b = np.array((
        2,
        7,
        3,
    ))

    c = np.array((
        -1,
        -2,
         0,
         0,
         0,
         0,
    ))

    init = np.array((2, 3, 4))

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