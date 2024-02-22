import numpy as np

def newton_raphson(f, x0, jacobian, tol=1e-6, max_iter=100, suppress_output=False):
    xk = x0
    iter_results = []

    for i in range(max_iter):
        fk = np.array([fi(*xk) for fi in f])
        jk = np.array([f_ij(*xk) for f_ij in jacobian])
        error = np.linalg.norm(fk)

        if not suppress_output:
            print('{:2} | {} | {} | {:.8f}'.format(i, xk, fk, error))

        xk -= np.linalg.solve(jk.reshape((fk.shape[0]), -1), fk)

        iter_results.append((xk, fk, error))

        if error <= tol:
            if not suppress_output:
                print('Root identified at x = {}'.format(xk))
            return iter_results

    if not suppress_output:
        print('No root identified!')
    return iter_results


if __name__ == '__main__':
    np.set_printoptions(precision=4, suppress=True, sign=' ', floatmode='fixed')

    # Scalar Test
    #
    F = [
        lambda x: x**4 - 7*x**3 + 17*x**2 - 17*x + 6,
    ]
    dF = [
        lambda x: 4*x**3 - 21*x**2 + 34*x - 17,
    ]
    res = newton_raphson(F, (1.1,), dF, tol=1e-10)
    assert(len(res))
    assert(np.isclose(res[-1][0], 1))

    res = newton_raphson(F, (2.2,), dF)
    assert(len(res))
    assert(np.isclose(res[-1][0], 2))

    # Multivariate Test
    #
    F = [
        lambda x1, x2: 3*x1*x2 + 7*x1 + 2*x2 - 3,
        lambda x1, x2: 5*x1*x2 - 9*x1 - 4*x2 + 6,
    ]
    dF = [
        lambda x1, x2: 3*x2 + 7,
        lambda x1, x2: 3*x1 + 2,
        lambda x1, x2: 5*x2 - 9,
        lambda x1, x2: 5*x1 - 4
    ]
    res = newton_raphson(F, (1, 2), dF, suppress_output=True)
    assert(len(res))
    assert(np.isclose(res[-1][0], np.array((0, 1.5))).all())
    assert(np.isclose(res[-1][1], np.array((0, 0))).all())

    print('All tests passed!')