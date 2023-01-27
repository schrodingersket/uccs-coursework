g1 = lambda x: 1 - x[0]**2 - x[1]**2
g2 = lambda x: 2**(1/2) - x[0] - x[1]
g3 = lambda x: x[1]

xa = [1/2, 1/2];
print('a: g1={:.2f} g2={:.2f} g3={:.2f}'.format(g1(xa), g2(xa), g3(xa)))

xb = [1, 0];
print('b: g1={:.2f} g2={:.2f} g3={:.2f}'.format(g1(xb), g2(xb), g3(xb)))

xc = [-1, 0];
print('c: g1={:.2f} g2={:.2f} g3={:.2f}'.format(g1(xc), g2(xc), g3(xc)))

xd = [-1/2, 0];
print('d: g1={:.2f} g2={:.2f} g3={:.2f}'.format(g1(xd), g2(xd), g3(xd)))

xe = [(1/2)**(1/2), (1/2)**(1/2)];
print('e: g1={:.2f} g2={:.2f} g3={:.2f}'.format(g1(xe), g2(xe), g3(xe)))