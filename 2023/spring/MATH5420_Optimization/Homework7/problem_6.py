import numpy as np

from quasi_newton import quadratic_sr1

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

res = quadratic_sr1(
    Q, 
    c, 
    np.array((0., 0., 0.)),
    B0,
    max_iter=216, 
    tol=10e-10, 
    suppress_output=False,
)