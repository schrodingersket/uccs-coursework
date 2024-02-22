import numpy as np

df = [
  lambda x1, x2, x3: 2*x1 + 2*x1*x3**2 + 2*x2,
  lambda x1, x2, x3: 2*x1 + 4*x2**3 + 8,
  lambda x1, x2, x3: 2*x1**2*x3,
]

Z = np.array((
  (-2, -1),
  ( 1,  0),
  (-1,  2),
))

print('(i)')
dfi = np.array([d(0, 0, 2) for d in df])
print('df = {}'.format(dfi))
print('Z df = {}'.format(np.matmul(Z.T, dfi)))
print('')

print('(ii)')
dfii = np.array([d(0, 0, 3) for d in df])
print('df = {}'.format(dfii))
print('Z df = {}'.format(np.matmul(Z.T, dfii)))
print('')

print('(iii)')
dfiii = np.array([d(1, 0, 1) for d in df])
print('df = {}'.format(dfiii))
print('Z df = {}'.format(np.matmul(Z.T, dfiii)))
print('')
