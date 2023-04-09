import numpy as np

def quadratic_gradient_descent(Q, c, x0, tol=1e-10, max_iter=10**5, suppress_output=False, x_star=None):
    """
    Implements gradient descent for a quadratic function of the form f(x) = .5 x^T Q x - c^T x.
    """
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

        gradient_norm = np.linalg.norm(jk)

        if gradient_norm < tol:
            if not suppress_output:
                print('Solution f(x) = {} identified at x = {}'.format(f(xk), xk))
            return iter_results

        # Step size (exact line search for quadratic functions)
        #
        step_numerator = np.dot(pk, jk)

        if step_numerator:
            step_length = -step_numerator / np.dot(pk, np.matmul(Q, pk))
        else:
            step_length = 0.


        # Newton step (x_{k+1})
        #
        xk_next = xk + step_length * pk

        # Compute observed rate constant when solution is known
        #
        if x_star is not None:
            fx_star = f(x_star)
            rate_constant_obs = np.abs(np.divide(f(xk_next) - fx_star, fk - fx_star))
        else:
            rate_constant_obs = None

        if not suppress_output:
            print('{:2} | {} | {:.6f} | {:.6f} | {:.6f}'.format(
                i, 
                xk, 
                fk, 
                gradient_norm, 
                rate_constant_obs or 0.,
            ))

        xk = xk_next

        iter_results.append((xk_next, fk, gradient_norm, rate_constant_obs))

    if not suppress_output:
        print('No root identified!')

    return iter_results


def quadratic_sr1(Q, c, x0, B0, tol=1e-10, max_iter=10**5, suppress_output=False, x_star=None):
    """
    Implements symmetric rank-one quasi-Newton method for a quadratic function of the form f(x) = .5 x^T Q x - c^T x.
    """
    xk = x0
    Bk = B0

    # Quadradic function and its gradient
    #
    f = lambda x: 0.5 * np.dot(x, np.matmul(Q, x)) - np.dot(c, x)
    jacobian = lambda x: np.matmul(Q, x) - c

    iter_results = [(x0, f(x0), np.linalg.norm(jacobian(x0)), jacobian(x0), None, None, B0)]

    for i in range(max_iter):
        # Iterate value
        #
        fk = f(xk)

        # Jacobian
        #
        jk = jacobian(xk)
        pk = np.linalg.solve(Bk, -jk)

        gradient_norm = np.linalg.norm(jk)

        if gradient_norm < tol:
            if not suppress_output:
                print('Solution f(x) = {} identified at x = {}'.format(f(xk), xk))
            return iter_results

        # Step size (exact line search for quadratic functions)
        #
        step_numerator = np.dot(pk, jk)

        if step_numerator:
            step_length = -step_numerator / np.dot(pk, np.matmul(Q, pk))
        else:
            step_length = 0.

        
        # Newton step (x_{k+1})
        #
        xk_next = xk + step_length * pk

        if not suppress_output:
            print('{:2} | {} | {:.6f} | {:.6f}'.format(
                i, 
                xk, 
                fk, 
                gradient_norm, 
            ))

        # Compute B_{k+1}
        #
        sk = xk_next - xk
        yk = jacobian(xk_next) - jk

        r1 = yk - np.matmul(Bk, sk)
        Bk = Bk + np.outer(r1, r1.T) / np.dot(r1.T, sk)

        xk = xk_next

        iter_results.append((xk_next, fk, gradient_norm, jk, sk, yk, Bk))

    if not suppress_output:
        print('No root identified!')

    return iter_results


def quadratic_bfgs(Q, c, x0, B0, tol=1e-10, max_iter=10**5, suppress_output=False, x_star=None):
    """
    Implements BFGS quasi-Newton method for a quadratic function of the form f(x) = .5 x^T Q x - c^T x.
    """
    xk = x0
    Bk = B0

    # Quadradic function and its gradient
    #
    f = lambda x: 0.5 * np.dot(x, np.matmul(Q, x)) - np.dot(c, x)
    jacobian = lambda x: np.matmul(Q, x) - c

    iter_results = [(x0, f(x0), np.linalg.norm(jacobian(x0)), jacobian(x0), None, None, B0)]

    for i in range(max_iter):
        # Iterate value
        #
        fk = f(xk)

        # Jacobian
        #
        jk = jacobian(xk)
        pk = np.linalg.solve(Bk, -jk)

        gradient_norm = np.linalg.norm(jk)

        if gradient_norm < tol:
            if not suppress_output:
                print('Solution f(x) = {} identified at x = {}'.format(f(xk), xk))
            return iter_results

        # Step size (exact line search for quadratic functions)
        #
        step_numerator = np.dot(pk, jk)

        if step_numerator:
            step_length = -step_numerator / np.dot(pk, np.matmul(Q, pk))
        else:
            step_length = 0.

        
        # Newton step (x_{k+1})
        #
        xk_next = xk + step_length * pk

        if not suppress_output:
            print('{:2} | {} | {:.6f} | {:.6f}'.format(
                i, 
                xk, 
                fk, 
                gradient_norm, 
            ))

        # Compute B_{k+1}
        #
        sk = xk_next - xk
        yk = jacobian(xk_next) - jk

        Bk = Bk - np.outer(np.matmul(Bk, sk), np.matmul(Bk, sk).T) / np.dot(sk.T, np.matmul(Bk, sk)) + np.outer(yk, yk) / np.dot(yk.T, sk)

        xk = xk_next

        iter_results.append((xk_next, fk, gradient_norm, jk, sk, yk, Bk))

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

    res = quadratic_gradient_descent(Q, c, (0, 0, 0), max_iter=216, tol=10e-8, suppress_output=True)

    assert(len(res))
    assert(np.isclose(res[-1][0], np.array((-1, -1/5, -1/25))).all())

    # Griva, Nash Sofer Example 12.9 (The Symmetric Rank-One Formula)
    #
    Q = np.array((
        (2, 0, 0),
        (0, 3, 0),
        (0, 0, 4),
    ))
    c = np.array((
        -8,
        -9,
        -8,
    ))

    res = quadratic_sr1(Q, c, np.array((0., 0., 0.)), np.eye(3), max_iter=216, tol=10e-8, suppress_output=False)

    assert(len(res) == 4)
    assert(np.isclose(res[-1][0], np.array((-4, -3, -2))).all())

    # Griva, Nash Sofer Example 12.11 (BFGS Formula)
    #
    Q = np.array((
        (2, 0, 0),
        (0, 3, 0),
        (0, 0, 4),
    ))
    c = np.array((
        -8,
        -9,
        -8,
    ))

    res = quadratic_bfgs(Q, c, np.array((0., 0., 0.)), np.eye(3), max_iter=216, tol=10e-8, suppress_output=False)

    assert(len(res) == 4)
    assert(np.isclose(res[-1][0], np.array((-4, -3, -2))).all())

    print('All tests passed!')