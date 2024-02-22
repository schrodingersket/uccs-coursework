import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

plt.title('Feasible Region')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

x_min, x_max = (-2, 50)
y_min, y_max = (-2, 50)

sample_size = 300

xx1, xx2 = np.meshgrid(
    np.linspace(x_min, x_max, sample_size),
    np.linspace(y_min, y_max, sample_size),
)

x1_min = 0
x2_min = 0

a0 = np.ones(xx1.shape) * x1_min
a1 = np.ones(xx1.shape) * x2_min

a2 = (3/2) * xx1 + 15
a3 = 2*xx1 + 12

plt.imshow(
    (
            (xx1 >= a0)
            & (xx2 >= a1)
            & (xx2 <= a2)
            & (xx2 <= a3)
    ).astype(int),
    extent=(x_min, x_max, y_min, y_max),
    origin='lower',
    cmap='Greys',
    alpha=0.25
)

plt.plot(xx1[0, :], a2[0, :], label=r'$x_2 \le \frac{3}{2}x_1 + 15$')
plt.plot(xx1[0, :], a3[0, :], label=r'$x_2 \le 2x_1 + 12$')

plt.plot(xx1[0, :], a1[0, :], 'k--', label=r'$x_2 \ge {}$'.format(x2_min))
plt.axvline(x=x1_min, label=r'$x_1 \ge {}$'.format(x1_min), color='gray', linestyle='dashed')

objective_fn = lambda z, x1: z/7 + (5/7)*x1
for i in np.arange(y_min, 4*y_max, (y_max - y_min) / 2):
    if i == y_min:
        plt.plot(xx1[0, :], objective_fn(i, xx1[0, :]), 'k-.', alpha=0.1, label='Level Sets')
    else:
        plt.plot(xx1[0, :], objective_fn(i, xx1[0, :]), 'k-.', alpha=0.1)

plt.legend(loc='upper right')
plt.xlim((x_min, x_max))
plt.ylim((y_min, y_max))
plt.grid(True)

plt.savefig('problem_5i.png', bbox_inches='tight')
