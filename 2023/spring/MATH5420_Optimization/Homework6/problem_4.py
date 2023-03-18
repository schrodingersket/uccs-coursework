import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

f = lambda x1, x2: 5*x1**4 + 6*x2**4 - 6*x1**2 + 2*x1*x2 + 5*x2**2 + 15*x1 - 7*x2 + 13

F = [
    lambda x1, x2: 20*x1**3 - 12*x1 + 2*x2  + 15,
    lambda x1, x2: 24*x2**3 + 2*x1  + 10*x2 - 7
]

dF = [
    lambda x1, x2: 60*x1**2,
    lambda x1, x2: 2,
    lambda x1, x2: 2,
    lambda x1, x2: 72*x2**2 + 10
]

print('{:2} | {:25} | {:25} | {:21}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(71, '-'))

x, _, _ = newton_raphson(F, np.array((1., 1.)), dF)[-1]
x1, x2 = x

print('')
print('Minimum value at root is f = {}'.format(f(*x)))

print('Gradient at root is given by {}'.format([ f(x1, x2) for f in F ]))

print('Hessian at root is given by {}'.format([ df(x1, x2) for df in dF]))

