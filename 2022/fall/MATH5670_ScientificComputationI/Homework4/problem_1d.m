maxiter = 800;       % number of iterations to take
m = 39;
ax = 0;
bx = 1;
alpha = 0;
beta = 0;
f = @(x) ones(size(x));   % f(x) = 1

h = (bx - ax)/(m + 1);
omega = 2 / (1 + sin(pi*h));

% Compute errors for various matrix splitting methods
%
[sor_optimal_errs, sor_optimal_rhoG] = iter_bvp_Asplit('SOR', omega, maxiter, m, ax, bx, alpha, beta, f);
[sor_naive_errs, sor_naive_rhoG] = iter_bvp_Asplit('naive_SOR', omega, maxiter, m, ax, bx, alpha, beta, f);

% plot errors vs. iteration:
figure
axis([0 maxiter 10e-16 max(max(sor_optimal_errs), max(sor_naive_errs))])
title('Naive SOR Method Error','FontSize',15)
xlabel('Iteration')

hold on;
semilogy(sor_optimal_errs, 'DisplayName', 'SOR');
semilogy(sor_naive_errs, 'DisplayName', 'Naive SOR');
legend;
hold off;

print('-dpng', sprintf('problem_1d_naive_sor_matrix_splitting_error_%d_iterations', maxiter))