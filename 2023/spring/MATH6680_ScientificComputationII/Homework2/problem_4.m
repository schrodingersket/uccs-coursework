% p13.m - solve linear BVP u_xx = exp(4x), u(-1)=u(1)=0
N = 64;

[D,x] = cheb(N);
D2 = D^2;
D2 = D2(2:N,2:N); % boundary conditions
f = exp(4*x(2:N));
u = D2\f; % Poisson eq. solved here
u = [0;u;0];
clf

subplot(2, 1, 1)
plot(x,u,'.','markersize',16)

xx = -1:.01:1;
uu = bary(x, u, xx); % interpolate grid data
line(xx,uu,'linewidth',.8)

line(xx(2:end-1),uu(2:end-1), 'linewidth',.8)
grid on
exact = ( exp(4*xx) - sinh(4)*xx - cosh(4) )/16;
title(['Poisson Equation max err = ' num2str(norm(uu(2:end-1)-exact(2:end-1),inf))],'fontsize',12)

subplot(2, 1, 2)
ex51 = @(x) 1./(1+16*x.^2);
u = ex51(x);
uu = bary(x, u, xx); % interpolate grid data

error = norm(uu-ex51(xx),inf);

plot(x,u,'.','markersize',16)
line(xx(2:end-1),uu(2:end-1), 'linewidth',.8)

grid on
title(['Exercise 5.1 max err = ' num2str(norm(uu(2:end-1)-ex51(xx)(2:end-1),inf))],'fontsize',12)

print('problem_4', '-dpng')