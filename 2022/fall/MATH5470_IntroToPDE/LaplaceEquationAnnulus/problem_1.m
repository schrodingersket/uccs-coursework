r1 = 3;
r2 = 5;
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(r1,r2,16));

X = R.*cos(T);
Y = R.*sin(T);
rr = linspace(r1, r2);

ua = @(r) 4 * log(r/5) / log(3/5);
ub = @(r) 4 * log(r/3) / log(5/3);

figure;
Za = 4 * log(sqrt(X.^2+Y.^2)/5) / log(3/5);
surf(X,Y,Za);
title('u(3, $\theta$) = 4;\; u(5, $\theta$) = 0', 'interpreter', 'latex', 'VerticalAlignment', 'bottom');
colormap turbo;
colorbar();
view(2);
print('-dpng', 'problem_1a_heatmap');

figure;
plot(rr, ua(rr));
title('Radius vs. Temperature');
xlabel('Radius (r)');
xlabel('Temperature (u)');

print('-dpng', 'problem_1a_u_vs_r');

figure;
Zb = 4 * log(sqrt(X.^2+Y.^2)/3) / log(5/3);
surf(X,Y,Zb);
title('u(3, $\theta$) = 0;\; u(5, $\theta$) = 4', 'interpreter', 'latex', 'VerticalAlignment', 'bottom');
colormap turbo;
colorbar();
view(2);
print('-dpng', 'problem_1b_heatmap');

figure;
plot(rr, ub(rr));
title('Radius vs. Temperature');
xlabel('Radius (r)');
xlabel('Temperature (u)');

print('-dpng', 'problem_1b_u_vs_r');

