import numpy as np

from plot_feasible_region import plot_2d_region

A = np.array((
    (1, 1, 1),
))

b = np.array((
    2,
))

c = np.array((
    1,
   -1,
    0,
))


plot_2d_region(
    A,
    b,
    c[0:2],
    xlims=(-2, 5),
    ylims=(-2, 5),
    usetex=True,
    legend=True,
    savefile='problem_5a.png',
    latex_vars=['x_1\'', 'x_2\''],
    plaintext_vars=['x1\'', 'x2\'']
)
