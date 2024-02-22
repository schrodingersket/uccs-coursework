import numpy as np

ddf = [
  lambda x1, x2, x3: 2 + 2*x3**2,
  lambda x1, x2, x3: 2,
  lambda x1, x2, x3: 4*x1*x3,

  lambda x1, x2, x3: 2,
  lambda x1, x2, x3: 12*x2**2,
  lambda x1, x2, x3: 0,

  lambda x1, x2, x3: 4*x1*x3,
  lambda x1, x2, x3: 0,
  lambda x1, x2, x3: 2*x1**2,
]

Z = np.array((
  (-2, -1),
  ( 1,  0),
  (-1,  2),
))

hessian = np.array([d(1, 0, 1) for d in ddf]).reshape(3, 3)
print('Hessian: {}'.format(hessian))

reduced_hessian = np.matmul(Z.T, np.matmul(hessian, Z))
print('Reduced Hessian: {}'.format(reduced_hessian))

eig, _ = np.linalg.eig(reduced_hessian)
print('Eigenvalues: {}'.format(eig))
