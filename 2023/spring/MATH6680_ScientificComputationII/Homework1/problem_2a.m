u = @(x) exp(sin(x).^2);
uprime = @(x) 2 * sin(x) .* cos(x) .* u(x);

clf;

[Nvec_fd, Nerror_fd, T_fd] = fd_diff_matrix(u, uprime);
subplot('position', [.1 .4 .8 .5]);

subplot(2, 1, 1);
hold on;
grid on;

loglog(Nvec_fd, Nerror_fd, '.', 'markersize', 15, 'DisplayName', 'Finite Difference Error')
semilogy(Nvec_fd, Nvec_fd.^(-4), '--', 'DisplayName', 'N^{-4}')

xlabel('N')
ylabel('Error')
title('Convergence of finite differences')
legend('Location', 'southwest');

[Nvec_spec, Nerror_spec, T_fd] = spectral_diff_matrix(u, uprime);

subplot(2, 1, 2);
hold on;
grid on;

loglog(Nvec_spec, Nerror_spec, '.', 'markersize', 15, 'DisplayName', 'Spectral Differentiation Error')

xlabel('N')
ylabel('Error')
title('Convergence of spectral differentiation')
legend('Location', 'southwest');

print('-dpng', 'problem_2a')