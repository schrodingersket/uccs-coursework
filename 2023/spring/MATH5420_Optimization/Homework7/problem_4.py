import numpy as np

from gradient_descent import quadratic_descent

np.set_printoptions(precision=6, suppress=True, sign=' ', floatmode='fixed')

Q = np.array((
    (8, 4),
    (4, 4),
))

c = np.array((
    3,
    0
))

res = quadratic_descent(
    Q, 
    c, 
    np.array((2., 2.)),
    max_iter=216, 
    tol=10e-10, 
    suppress_output=False,
)