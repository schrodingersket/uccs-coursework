import numpy as np

f = lambda x1, x2: 3 * (x1**4) - 2*(x1**3)*x2 - 4*(x1**2)*(x2**2) + 5*x1*(x2**3) + 2*(x2**4)

grad_f = (
  lambda x1, x2: 12*(x1**3) - 6*(x1**2)*x2 - 8*x1*(x2**2) + 5*(x2**3),
  lambda x1, x2: -2*(x1**3) - 8*(x1**2)*x2 + 15*x1*(x2**2) + 8*(x2**3)
)

hessian_f = (
  lambda x1, x2: 36*(x1**2) - 12*x1*x2 - 8*(x2**2),
  lambda x1, x2: -6*(x1**2) - 16*x1*x2 + 15*(x2**2),
  lambda x1, x2: -6*(x1**2) - 16*x1*x2 + 15*(x2**2),
  lambda x1, x2: -8*(x1**2) + 30*x1*x2 + 24*(x2**2)
)

x = np.array([1, -1])
p = np.array([0.1, 0.01])
print('x = {}'.format(x))
print('p = {}'.format(p))

# f(x)
#
print('f(x): {}'.format(f(*x)))

grad_f0 = np.array(
  [grad_f[0](*x), grad_f[1](*x)],
)
print('D f(x): {}'.format(grad_f0))

hessian_f0 = np.array(
  [
    [hessian_f[0](*x), hessian_f[1](*x)],
    [hessian_f[2](*x), hessian_f[3](*x)],
  ]
)

print('D2 f(x): {}'.format(hessian_f0.flatten()))

print('\n')

# f(x + p)
#
print('f(x + p): {}'.format(f(*(x+p))))

taylor_approx = f(*x) + np.matmul(p.T, grad_f0) + .5 * np.matmul(p.T, np.matmul(hessian_f0, p))
print('F(x + p): {} + {} + {} = {}'.format(
  f(*x), np.matmul(p.T, grad_f0), .5 * np.matmul(p.T, np.matmul(hessian_f0, p)),
  taylor_approx
))

print('|F(x+p) - f(x+p)| = {}'.format(np.abs(taylor_approx - f(*(x+p)))))
print('|p|^2: {}'.format(np.linalg.norm(p, 1) ** 2))