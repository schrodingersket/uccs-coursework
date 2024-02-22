import numpy as np

df_x0 = np.array((
  4,
  4
))

p0 = -1 * df_x0

Q = np.array((
  (6, -2),
  (-2, 6)
))

alpha0 = -np.dot(df_x0, p0) / np.dot(p0.T, np.matmul(Q, p0))
print('alpha0 = -1 * {}/{} = {}'.format(np.dot(df_x0, p0), np.dot(p0.T, np.matmul(Q, p0)), alpha0))
print('eig(Q) = {}'.format(np.linalg.eig(Q)))

