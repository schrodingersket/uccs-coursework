import numpy as np

def quadratic_descent(Q, c, x0, tol=1e-10, max_iter=10**5, suppress_output=False):
    xk = x0
    iter_results = []

    # Quadradic function and its gradient
    #
    f = lambda x: 0.5 * np.dot(x, np.matmul(Q, x)) - np.dot(c, x)
    jacobian = lambda x: np.matmul(Q, x) - c
    fk = f

    for i in range(max_iter):
        # Iterate value
        #
        fk = f(xk)

        # Jacobian
        #
        jk = jacobian(xk)
        pk = -jk

        # Step size (exact line search for quadratic functions)
        #
        step_length = -np.dot(pk, jk) / np.dot(pk, np.matmul(Q, pk))

        error = np.linalg.norm(jk)

        if not suppress_output:
            print('{:2} | {} | {} | {:.8f}'.format(i, xk, fk, error))

        xk += step_length * pk

        iter_results.append((xk, fk, error))

        if error < tol:
            if not suppress_output:
                print('Solution identified at x = {}'.format(xk))
            return iter_results

    if not suppress_output:
        print('No root identified!')

    return iter_results


if __name__ == '__main__':
    np.set_printoptions(precision=4, suppress=True, sign=' ', floatmode='fixed')

    # Griva, Nash, Sofer Example 12.1
    #
    Q = np.array((
        (1, 0, 0),
        (0, 5, 0),
        (0, 0, 25),
    ))
    c = np.array((
        -1,
        -1,
        -1,
    ))


    res = quadratic_descent(Q, c, (0, 0, 0), max_iter=216, tol=10e-8, suppress_output=False)

    assert(len(res))
    assert(np.isclose(res[-1][0], np.array((-1, -1/5, -1/25))).all())

    print('All tests passed!')