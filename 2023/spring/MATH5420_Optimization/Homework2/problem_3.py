import numpy as np

from newton_raphson import newton_raphson

np.set_printoptions(precision=4, suppress=True, formatter={
    'float': lambda x: '{0:0.8f}'.format(x),
})

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

newton_raphson(F, (1., 2.), dF)