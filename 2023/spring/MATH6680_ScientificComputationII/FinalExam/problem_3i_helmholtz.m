clf;
clear all;

% Gaussian RBF
%
rbf = @(e, r) exp(-(e*r).^2);
d2rbf = @(e, r) 2 * e^2 * (2 * (e*r).^2 - 1) .* rbf(e, r);

exact_eigenvalues = @(j) (pi^2 * j.^2) ./ 4;

% Eigenvalues to compute
%
M = 5;

% Collocation points
%
N = 24;

[D2, x] = D2RBF(N, rbf, d2rbf);
D2 = D2(2:N, 2:N);

[V, lam] = eig(-D2);
[eigenvalues, ii] = sort(diag(lam), 'ascend');
wave_numbers = sqrt(eigenvalues);
ii = ii(1:M);
eigenvectors = real(V(:, ii));

for k=1:M
    fprintf('lambda_%d: \t %14.11f \t %14.11f \n', k, eigenvalues(k), norm(eigenvalues(k) - exact_eigenvalues(k)))
end

% Plot eigenvalues
%
subplot(2, 3, 1)
semilogy(1:M, eigenvalues(1:M), '.', 'markersize', 12)
grid on
axis square
title(['N = ' int2str(N) ', \lambda_{min} = ' num2str(min(real(eigenvalues(1:M))), '%15.3f')]);

% Plot eigenvectors
%
xx = linspace(-1, 1);
for i = 1:M
    subplot(2, 3, i+1)
    plot(xx, interp1(x, [0; eigenvectors(:, i); 0], xx), 'linewidth', 2)
    grid on
    axis square
    title(['\lambda = ' num2str(eigenvalues(i), '%15.3f')]);
end

print('-dpng', 'problem_3i_helmholtz.png')