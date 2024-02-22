import numpy as np
import matplotlib.pyplot as plt

from basic_solutions import compute_basic_solutions

# Plot constraints
#
plt.rcParams['text.usetex'] = True

plt.title('Feasible Region')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

x_min, x_max = (-.5, 6)
y_min, y_max = (-.5, 6)

sample_size = 300

xx1, xx2 = np.meshgrid(
    np.linspace(x_min, x_max, sample_size),
    np.linspace(y_min, y_max, sample_size),
)

x1_min = 0
x2_min = 0

a0 = np.ones(xx1.shape) * x1_min
a1 = np.ones(xx1.shape) * x2_min

a2 = 4 - (4 * xx1) / 3
a3 = (1/2) * xx1 - 1

plt.imshow(
    (
            (xx1 >= a0)
            & (xx2 >= a1)
            & (xx2 <= a2)
            & (xx2 >= a3)
    ).astype(int),
    extent=(x_min, x_max, y_min, y_max),
    origin='lower',
    cmap='Greys',
    alpha=0.25
)

plt.plot(xx1[0, :], a2[0, :], label=r'$x_2 \le -\frac{4}{3}x_1 + 4$')
plt.plot(xx1[0, :], a3[0, :], label=r'$x_2 \ge \frac{1}{2}x_1 - 1$')

plt.plot(xx1[0, :], a1[0, :], 'k--', label=r'$x_2 \ge {}$'.format(x2_min))
plt.axvline(x=x1_min, label=r'$x_1 \ge {}$'.format(x1_min), color='gray', linestyle='dashed')

plt.plot([1], [1], marker='o', markersize=3, markerfacecolor='green')

plt.legend(loc='upper right')
plt.xlim((x_min, x_max))
plt.ylim((y_min, y_max))
plt.grid(True)

plt.savefig('problem_9.png', bbox_inches='tight')


# Find basic solutions
#
A = np.array((
    (4, 3, 1, 0),
    (1, -2, 0, 1),
))

b = np.array((
    12,
    2,
))

rows, cols = A.shape
solns = compute_basic_solutions(A, b, var_names=['x1', 'x2', 's1', 's2'])

for s in solns:
    soln = s.get('solution')
    pretty_basis = s.get('formatted_basis')
    print('Basis: {}'.format(pretty_basis))
    print('B: {}'.format(np.array([A[:, k] for k in s.get('basis')]).T))

    if soln is not None:
        full_solution = [0] * cols
        for i, beta in enumerate(s.get('basis')):
            full_solution[beta] = soln[i]

        feasible = all(i >= 0 for i in full_solution)
        print('Solution ({}feasible): {}'.format('' if feasible else 'in', full_solution))
    else:
        print('No solution exists!')

    print('')
    input('Press [Enter] to see the next point:')