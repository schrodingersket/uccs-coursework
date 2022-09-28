maxiter = 100;       % number of iterations to take
m = 39;
ax = 0;
bx = 1;
alpha = 0;
beta = 0;
f = @(x) ones(size(x));   % f(x) = 1

h = (bx - ax)/(m + 1);
omega = linspace(1e-16, 2);
rhoG = zeros(size(omega));

% Compute spectral radius for various values of omega
%
for i=1:length(omega)
[_, rhoG(i)] = iter_bvp_Asplit('SOR', omega(i), maxiter, m, ax, bx, alpha, beta, f);
end


% plot errors vs. iteration:
figure
axis([0 2 min(rhoG) max(rhoG)])
title('Iteration Matrix Spectral Radii', 'FontSize',15)
xlabel('\omega')
ylabel('\rho(G(\omega))')

hold on;
plot(omega, rhoG);
hold off;

print('-dpng', 'problem_1c_spectral_radius')