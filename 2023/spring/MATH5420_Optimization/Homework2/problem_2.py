import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

F = [
    lambda x: x**3 - 5*x**2 - 12*x + 19,
]

dF = [
    lambda x: 3*x**2 - 10*x - 12,
]

print('{:2} | {:15} | {:15} | {:13}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(51, '-'))

newton_raphson(F, np.array((-3.,)).reshape(1,1), dF)
print('')

newton_raphson(F, np.array((1.,)).reshape(1,1), dF)
print('')

newton_raphson(F, np.array((8.,)).reshape(1,1), dF)