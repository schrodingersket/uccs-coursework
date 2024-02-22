import numpy as np

c = np.array((
  -1,
   3
))

Q = np.array((
  (2, -1),
  (-1, 3)
))

print('x_* = {}'.format(-1/np.sqrt(np.dot(c.T, (np.matmul(np.linalg.inv(Q), c)))) * np.matmul(np.linalg.inv(Q), c)))

print('-1/sqrt(3) = {}'.format(-1/(3**(0.5))))
