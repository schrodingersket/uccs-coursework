import numpy as np

def newton_raphson(f, x0, jacobian, tol=1e-6, max_iter=100, suppress_output=False):
    """
    Returns an array of (xk, fk, error)
    """
    xk = x0
    iter_results = []

    for i in range(max_iter):
        fk = np.array([fi(*xk) for fi in f])
        jk = np.array([f_ij(*xk) for f_ij in jacobian])
        error = np.linalg.norm(fk)

        if not suppress_output:
            print('{:2} | {} | {} | {:.8f}'.format(i, xk, fk, error))

        newton_step = -1 * np.linalg.solve(jk.reshape((fk.shape[0]), -1), fk)
        xk += newton_step

        iter_results.append((xk, fk, error, newton_step))

        if error <= tol:
            if not suppress_output:
                print('Root identified at x = {}'.format(xk))
            return iter_results

    if not suppress_output:
        print('No root identified!')
    return iter_results


if __name__ == '__main__':
    np.set_printoptions(precision=4, suppress=True, sign=' ', floatmode='fixed')

    # Griva, Nash Sofer Example 11.5
    #
    F = [
        lambda x1, x2: 10*x1 - 3*x2,
        lambda x1, x2: 14*x2 - 3*x1,
    ]
    dF = [
        lambda x1, x2: 10,
        lambda x1, x2: -3,
        lambda x1, x2: -3,
        lambda x1, x2: 14,
    ]
    res = newton_raphson(F, np.array((2., 3.)), dF, suppress_output=False)
    assert(len(res))

    print('All tests passed!')