import numpy as np
import array_to_latex as a2l

from quasi_newton import quadratic_bfgs

to_ltx = lambda arr: a2l.to_ltx(arr, frmt = '{:6.5f}', arraytype = 'bmatrix', print_out=False)

np.set_printoptions(precision=6, suppress=True, sign=' ', floatmode='fixed')

Q = np.array((
    (5, 2, 1),
    (2, 7, 3),
    (1, 3, 9),
))

c = np.array((
   -9,
    0,
   -8
))

B0 = np.eye(3)

res = quadratic_bfgs(
    Q, 
    c, 
    np.array((0., 0., 0.)),
    B0,
    max_iter=216, 
    tol=10e-10, 
    suppress_output=False,
)

for i, (xk_next, fk, gradient_norm, jk, jk_next, sk, yk, Bk, pk, step_length) in enumerate(res):
    print('Iteration {}:'.format(i))
    print('div(f)_{}: {}'.format(i, to_ltx(jk)))
    print('||div(f)_{}||: {}'.format(i, np.linalg.norm(jk)))
    print('p_{}: {}'.format(i, to_ltx(pk)))
    print('alpha_{}: {}'.format(i, step_length))

    print('')
    print('x_{}: {}'.format(i+1, to_ltx(xk_next)))
    print('div(f)_{}: {}'.format(i+1, to_ltx(jk_next)))
    print('s_{}: {}'.format(i, to_ltx(sk)))
    print('y_{}: {}'.format(i, to_ltx(yk)))

    print('Bk_{}: {}'.format(i+1, to_ltx(Bk)))
    print('')
    print('')

