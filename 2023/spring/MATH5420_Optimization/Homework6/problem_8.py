import numpy as np
import scipy.linalg

A = np.array((
    (1, 4, 3),
    (4, 2, 5),
    (3, 5, 3),
))

E = 8 * np.identity(3)

lu, d, perm = scipy.linalg.ldl(A + E)

print(lu)
print(d)
print(perm)