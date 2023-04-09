import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

# Griva, Nash Sofer Example 11.5
#
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
res = newton_raphson(F, x0, dF, suppress_output=True)

(xk, dfk, error, pk) = res[0]

alpha = 1.
mu = 0.2

f = lambda x: np.array([fi(*x) for fi in F])
df = lambda x: np.array([f_ij(*x) for f_ij in dF])

ff = f(x0 + alpha * pk)
print(ff)

f_armijo = f(x0) + mu * alpha * np.matmul(pk.T, df(x0).reshape((ff.shape[0]), -1).T)
print(f_armijo)
