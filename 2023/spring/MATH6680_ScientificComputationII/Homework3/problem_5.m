% problem_5.m - eignmodes of u_xxxx + u_xxx = lambda u_xx on \omega=[-2, 2] with
%               u(+/- 2) = u_x(+/- 2) = 0
%
clear all
close all

N = 50;
M = 5; % number of eigenvalues to compute
[D, x] = cheb(N);


D2 = (D^2)(2:N, 2:N);

S = diag([0; 1 ./ (1 - x(2:N).^2); 0]);

D3 = ((diag(1 - x.^2) * D^3 - 6 * diag(x) * D^2 - 6 * D)*S)(2:N, 2:N);
D4 = ((diag(1 - x.^2) * D^4 - 8 * diag(x) * D^3 - 12 * D^2)*S)(2:N, 2:N);

% With our variable transformation x = 2t, our system (in t) becomes:
%
%  (1/2)^4 u_tttt + (1/2)^3 u_ttt = lambda (1/2)^2 u_tt on \omega'=[-1, 1] with
%   u(+/- 1) = u_t(+/- 1) = 0
%
A = (1/2)^4 * D4 + (1/2)^3 * D3;
B = (1/2)^2 * D2;

[V, lam] = eig(A, B);
[eigenvalues, ii] = sort(diag(lam), 'descend');
ii = ii(1:M);
eigenvectors = real(V(:, ii));

% Plot eigenvalues
%
semilogy(1:M, eigenvalues(1:M), '.', 'markersize', 12)
grid on

axis square
title(['N = ' int2str(N) ', \lambda_{max} = ' num2str(max(real(eigenvalues(1:M))), '%15.11f')]);
drawnow

print('-dpng', 'problem_5_eigenvalues.png')

% Plot eigenvectors
%
figure
disp(eigenvalues(1:M))
xx = linspace(-1, 1);
for i = 1:M
    subplot(2, 3, i)
    plot(xx, interp1(x, [0; eigenvectors(:, i); 0], xx), 'linewidth', 2)
    grid on
    axis square
    title(['\lambda = ' num2str(eigenvalues(i), '%15.11f')]);
end

print('-dpng', 'problem_5_eigenvectors.png')