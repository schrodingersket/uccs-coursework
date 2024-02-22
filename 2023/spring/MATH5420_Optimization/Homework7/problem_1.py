import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

F = [
    lambda x1, x2: 2*(x1 - 2*x2) + 4*(x1**3),
    lambda x1, x2: 2*(x1 - 2*x2)*(-2),
]
dF = [
    lambda x1, x2: 2 + 12*(x1**2),
    lambda x1, x2: -4,
    lambda x1, x2: -4,
    lambda x1, x2: 8,
]
x0 = np.array((2., 1.))

alpha = 1.
mu = 0.2

f = lambda x1, x2: (x1 - 2*x2)**2 + x1**4
df = lambda x: np.array([fi(*x) for fi in F])
ddf = lambda x: np.array([f_ij(*x) for f_ij in dF])

gradient = df(x0)
hessian = ddf(x0).reshape(x0.shape[0], -1)
p = np.linalg.solve(-hessian, gradient)  # Newton step

print('x + alpha p = {}'.format(x0 + alpha * p))
print('p = {}'.format(p))
print('f(x0 + alpha * p) = {}'.format(f(*(x0 + alpha * p))))
print('f(x0) + mu * alpha * p^T div(f(x0) = {}'.format(f(*x0) + mu * alpha * np.dot(p.T, gradient)))
print('f(x0) = {}'.format(f(*x0)))
print('div(f(x0)) = {}'.format(df(x0)))

# Ensure we computed gradient and Hessian correctly by checking Newton's method convergence
#
res = newton_raphson(F, x0, dF, suppress_output=True)

assert(len(res) < 100)