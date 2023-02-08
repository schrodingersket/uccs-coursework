import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

F = [
    lambda x1, x2: x1**2 + x2**2 - 1,
    lambda x1, x2: 5*x1**2 - x2 - 2,
]

dF = [
    lambda x1, x2: 2*x1,
    lambda x1, x2: 2*x2,
    lambda x1, x2: 10*x1,
    lambda x1, x2: -1
]

print('{:2} | {:23} | {:23} | {:21}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(71, '-'))

newton_raphson(F, np.array((1., 1.)), dF)
print('')

newton_raphson(F, np.array((-1., 1.)), dF)
print('')

newton_raphson(F, np.array((0.5, -0.5)), dF)
print('')

newton_raphson(F, np.array((-1., -1.)), dF)