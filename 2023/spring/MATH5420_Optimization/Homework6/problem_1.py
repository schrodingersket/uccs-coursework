import numpy as np

f = lambda x1, x2: 8*x1**2 + 3*x1*x2 + 7*x2**2 - 25*x1 + 31*x2 - 29


A = np.array((
    (16, 3),
    (3, 14),
))

b = np.array((
    25,
   -31
))

x = np.linalg.solve(A, b)

print('Stationary point at x = {}'.format(x))
print('Minimum of f = {}'.format(f(*x)))