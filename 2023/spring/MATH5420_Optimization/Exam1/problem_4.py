import numpy as np

f = lambda x: 0.5 * np.matmul(np.matmul(x.T, Q), x) - np.matmul(x, c.T)

Q = np.array((
    (1, 2, 2),
    (2, 4, -1),
    (2, -1, 4)
))

c = np.array((
    1,
    1,
    1
))

print('Eigenvalues of Q: {}'.format(np.linalg.eigvals(Q)))

print('inv(Q): {}'.format(np.linalg.inv(Q)))

# x_1 = np.matmul(np.linalg.inv(Q), c)
x_1 = np.linalg.solve(Q, c)

print('Single Newton iteration: {}'.format(x_1))

print('Gradient at x_1: {}'.format(np.matmul(Q, x_1) - c))

print('f(x_1): {}'.format(
    f(x_1)
))

x_min = x_1 + np.array((
    1.0,
   -1.0,
   -1.0
))
print('f({}): {}'.format(
    x_min,
    f(x_min)
))