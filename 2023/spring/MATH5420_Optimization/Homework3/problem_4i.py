import numpy as np

from basic_solutions import compute_basic_solutions

A = np.array((
    (2, 1, 1, 0),
    (1, 1, 0, 1)
))

b = np.array((
    100,
    80
))

compute_basic_solutions(A, b)
