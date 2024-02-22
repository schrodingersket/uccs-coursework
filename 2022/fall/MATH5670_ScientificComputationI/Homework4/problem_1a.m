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
[jacobi_errs, jacobi_rhoG] = iter_bvp_Asplit('Jacobi', omega, maxiter, m, ax, bx, alpha, beta, f);
[gs_errs, gs_rhoG] = iter_bvp_Asplit('GS', omega, maxiter, m, ax, bx, alpha, beta, f);
[sor_errs, sor_rhoG] = iter_bvp_Asplit('SOR', omega, maxiter, m, ax, bx, alpha, beta, f);

% plot errors vs. iteration:
figure
axis([0 maxiter 10e-16 1])
title('Matrix Splitting Method Errors','FontSize',15)
xlabel('Iteration')

hold on;
semilogy(jacobi_errs, 'DisplayName', 'Jacobi');
semilogy(gs_errs, 'DisplayName', 'Gauss-Seidel');
semilogy(sor_errs, 'DisplayName', 'SOR');
legend;
hold off;

print('-dpng', sprintf('problem_1a_matrix_splitting_error_%d_iterations', maxiter))