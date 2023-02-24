import numpy as np

from plot_feasible_region import plot_2d_region



plot_2d_region(
    np.array((
        (4, 1),
        (1, 1),
        (1, 0),
    )),
    np.array((
        100,
        80,
        40,
    )),
    xlims=(-2, 100),
    ylims=(-2, 100),
    usetex=True,
    legend=True,
    savefile='problem_1v.png',
    objective_fn=lambda z, x1: (1/8)*z - (7/8)*x1,  # Solve z = 7*x_1 + 8*x_2 for x_2
)