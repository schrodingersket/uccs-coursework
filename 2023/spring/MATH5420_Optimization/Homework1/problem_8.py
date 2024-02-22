import numpy as np
  
x0 = 0

def f(x):
  return np.log(1 + x)

def third_order_approx(x):
  return x - 1/2 * np.power(x, 2) + 1/3 * np.power(x, 3)

p1 = 0.1
p2 = 0.01

print('p = 0.1')
print('F(x + p) = {:.6f}'.format(third_order_approx(x0 + p1)))
print('f(x + p) = {:.6f}'.format(f(x0 + p1)))
print('|f(x + p) - F(x + p)| = {:.6f}'.format(
  np.abs(f(x0 + p1) - third_order_approx(x0 + p1)))
)
print(' ')

print('p = 0.01')
print('F(x + p) = {:.12f}'.format(third_order_approx(x0 + p2)))
print('f(x + p) = {:.12f}'.format(f(x0 + p2)))
print('|f(x + p) - F(x + p)| = {:.12f}'.format(
  np.abs(f(x0 + p2) - third_order_approx(x0 + p2)))
)