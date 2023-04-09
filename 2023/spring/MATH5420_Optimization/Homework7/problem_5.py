import numpy as np

from quasi_newton import quadratic_gradient_descent

np.set_printoptions(precision=6, suppress=True, sign=' ', floatmode='fixed')

gammas = [1, 10, 100, 1000]

x_stars = [
    np.array((1., 1., 1.)),
    np.array((1., 0.1, 0.01)),
    np.array((1., 0.01, 0.0001)),
    np.array((1., 0.001, 0.000001)),
]

Q = lambda gamma: np.array((
    (1,     0,        0),
    (0, gamma,        0),
    (0,     0, gamma**2)
))


c = np.array((
    1,
    1,
    1
))

results = [quadratic_gradient_descent(
    Q(g), 
    c, 
    np.array((1.1, 0.5, 0.05)),
    max_iter=10**8, 
    tol=10e-10, 
    suppress_output=True,
    x_star=x_stars[i]
) for i, g in enumerate(gammas)]

for i, r in enumerate(results):
    cond = np.linalg.cond(Q(gammas[i]))
    rate_constant_bound = ((cond - 1)/(cond + 1)) ** 2
    print('Condition number: {}'.format(np.linalg.cond(Q(gammas[i]))))
    print('Rate constant bound: {}'.format(rate_constant_bound))

    print(' i |                 xk              |   f(xk)   |  error   | rate     | bound')

    enumerated_results = [(i, result) for i, result in enumerate(r)]
    [print('{:2} | {} | {:.6f} | {:.6f} | {:.6f} | {:.6f}'.format(
        i,
        xk,
        fk,
        gradient_norm,
        rate_constant or 0.,
        rate_constant_bound,
    )) for i, (xk, fk, gradient_norm, rate_constant) in enumerated_results[:10]]

    if len(r) > 10:
        print(' ...({} more iterations)'.format(len(r) - 10))
    print(' Solution f(x) = {} identified at x = {}'.format(r[-1][1], r[-1][0]))

    print('')