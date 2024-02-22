import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

f = lambda x1, : 5*x1**5 + 2*x1**3 - 4*x1**2 - 3*x1 + 2

F = [
    lambda x: 25*x**4 + 6*x**2 - 8*x - 3
]

dF = [
    lambda x: 100*x**3 + 12*x - 8
]

print('{:2} | {:25} | {:25} | {:21}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(71, '-'))

x, _, _ = newton_raphson(F, np.array((1.,)), dF)[-1]
x1, = x

print('Minimum value at root is f = {}'.format(f(*x)))

print('Gradient at root is given by {}'.format([ f(x1) for f in F ]))

print('Hessian at root is given by {}'.format([ df(x1) for df in dF]))