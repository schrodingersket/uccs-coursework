import os
import numpy as np
import matplotlib.pyplot as plt

x1 = np.linspace(0, 2, 100)

fig, ax = plt.subplots()

def g1(x):
    return np.sqrt(4 - np.power(x, 2))

def f(x):
    return x

# g1(x)
#
plt.plot(x1, g1(x1), color='silver', label=r'4 - $x_1^2 - x_2^2 \geq 0$')
plt.plot(x1, -1 * g1(x1), color='silver')

# g2(x)
#
plt.axvline(1, color='gray', label=r'$x_1^2 - 1 \geq 0$')

# f(x)
#
plt.plot(x1, f(x1), color='black', label=r'$f(x) = x_1$')

plt.xlim((-0.1, np.pi))
plt.ylim((-np.pi, np.pi))
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')

# Feasible region
#
plt.fill_between(
    x1,
    g1(x1),
    -g1(x1),
    where=x1 >= 1,
    color='green',
    alpha=0.25,
    label=r'$S$',
)

plt.legend()
plt.title(r'Feasible Region $S$ for $f(x)$')
plt.savefig('{}.jpg'.format(os.path.splitext(os.path.basename(__file__))[0]))