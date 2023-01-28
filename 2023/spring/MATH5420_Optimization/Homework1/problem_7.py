import numpy as np

w, v = np.linalg.eig(np.array([
    [4, -3],
    [-3, 10],
]))

print(w)