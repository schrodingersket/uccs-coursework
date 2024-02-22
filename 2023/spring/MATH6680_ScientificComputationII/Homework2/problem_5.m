% p14.m - solve nonlinear BVP u_xx + 4u_x + exp(x)u = sin(8x), u(-1)=u(1)=0
% (compare p13.m)
N = 16;
[D,x] = cheb(N);
D2 = D^2;
D2 = D2(2:N,2:N);
D4 = 4 * D(2:N,2:N);
u = zeros(N-1,1);
change = 1;
it = 0;
xN = x(2:N);
while change > 1e-11 % fixed-point iteration
  unew = (D2 + D4) \ (sin(8*xN) -  + exp(xN).*u);
  change = norm(unew-u,inf);
  u = unew;
  it = it+1;
  fprintf('u(0) = %.11f\n', unew(N/2))
end
u = [0;u;0];
clf
plot(x,u,'o','markersize', 4)
xx = -1:.01:1;
uu = polyval(polyfit(x,u,N),xx);
line(xx,uu,'linewidth',.8)
grid on
title('Solution to $u_{xx} + 4u_x + e^x u = \sin{(8x)}$', 'interpreter', 'latex')
xlabel('$x$', 'interpreter', 'latex')
xlabel('$u(x)$', 'interpreter', 'latex')

print('-dpng','problem_5')
