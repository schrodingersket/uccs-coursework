maxiter = 500;       % number of iterations to take
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
[gs_errs, gs_rhoG] = iter_bvp_Asplit('GS', omega, maxiter, m, ax, bx, alpha, beta, f);
[backward_gs_errs, backward_gs_rhoG] = iter_bvp_Asplit('backward_GS', omega, maxiter, m, ax, bx, alpha, beta, f);

% plot errors vs. iteration:
figure
axis([0 maxiter min(min(gs_errs), min(backward_gs_errs)) max(max(gs_errs), max(backward_gs_errs))])
title('Forward and backward Gauss-Seidel Error','FontSize',15)
xlabel('Iteration')

hold on;
semilogy(gs_errs, 'DisplayName', 'Forward Gauss-Seidel');
semilogy(backward_gs_errs, 'DisplayName', 'Backward Gauss-Seidel');
legend;
hold off;

print('-dpng', sprintf('problem_2b_forward_backward_gs_error_%d_iterations', maxiter))