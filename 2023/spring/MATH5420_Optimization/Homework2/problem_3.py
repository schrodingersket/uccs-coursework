import numpy as np
import matplotlib.pyplot as plt

from newton_raphson import newton_raphson

np.set_printoptions(precision=8, suppress=True, sign=' ', floatmode='fixed')

a = 2

print('a = {}:'.format(a))
print('')
print('{:2} | {:15} | {:15} | {:13}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(51, '-'))

iterations = newton_raphson([lambda x: x**2 - a], np.array((4.,)).reshape(1,1), [lambda x: 2*x])
print('')

fig, ax = plt.subplots(constrained_layout=True)
xx = np.linspace(0.3, len(iterations) - 1)
ax.plot(range(len(iterations)), [err for xk, fk, err in iterations], label='Actual Error')
ax.plot(xx, np.power(xx, -2), label='$x^{-2}$')
ax.legend()
ax.set_xlabel('Iterations')
ax.set_ylabel('Error')
ax.set_title('Newton-Raphson Method Error')
plt.savefig('problem_3.png', bbox_inches='tight')

a = -2

# Real initial guess
#
print('a = {}:'.format(a))
print('')
print('{:2} | {:15} | {:15} | {:13}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(51, '-'))
newton_raphson([lambda x: x**2 - a], np.array((1.,)).reshape(1,1), [lambda x: 2*x], max_iter=10)
print('')

# Complex initial guess
#
print('a = {}:'.format(a))
print('')
print('{:2} | {:27} | {:27} | {:13}'.format('i', 'xk', 'f(xk)', '||f(xk)||'))
print(''.ljust(75, '-'))

newton_raphson([lambda x: x**2 - a], np.array(1 + 1j).reshape(1,1), [lambda x: 2*x])