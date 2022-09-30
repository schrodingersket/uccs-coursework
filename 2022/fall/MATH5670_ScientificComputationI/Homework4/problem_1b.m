maxiter = 800;       % number of iterations to take
m = 39;
ax = 0;
bx = 1;
alpha = 0;
beta = 0;
f = @(x) ones(size(x));   % f(x) = 1

h = (bx - ax)/(m + 1);
omega = 2 / (1 + sin(pi*h));
epsilon = 0.05;

% Compute errors for various matrix splitting methods
%
[sor_optimal_errs, sor_optimal_rhoG] = iter_bvp_Asplit('SOR', omega, maxiter, m, ax, bx, alpha, beta, f);
[sor_high_errs, sor_high_rhoG] = iter_bvp_Asplit('SOR', omega + 2*epsilon, maxiter, m, ax, bx, alpha, beta, f);
[sor_low_errs, sor_low_rhoG] = iter_bvp_Asplit('SOR', omega - epsilon, maxiter, m, ax, bx, alpha, beta, f);

% plot errors vs. iteration:
figure
axis([0 maxiter 10e-16 1])
title('SOR Method Errors','FontSize',15)
xlabel('Iteration')

hold on;
semilogy(sor_optimal_errs, 'DisplayName', sprintf('\\omega = %0.2f', omega));
semilogy(sor_low_errs, 'DisplayName', sprintf('\\omega - \\epsilon = %0.2f', omega - epsilon));
semilogy(sor_high_errs, 'DisplayName', sprintf('\\omega + 2\\epsilon = %0.2f', omega + 2*epsilon));
legend;
hold off;

print('-dpng', sprintf('problem_1b_sor_matrix_splitting_error_%d_iterations', maxiter))