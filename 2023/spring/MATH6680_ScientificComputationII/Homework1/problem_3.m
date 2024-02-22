u = @(x) exp(sin(x));
uprime = @(x) cos(x) .* u(x);

clf;

[Nvec_fd, Nerror_fd, T_fd] = fd_diff_matrix(u, uprime);
subplot('position', [.1 .4 .8 .5]);

subplot(2, 1, 1);
hold on;
grid on;

h = (2*pi./Nvec_fd);
u5max = 24.811;
leading_taylor_error = (1/30) * h.^4 * u5max;
loglog(Nvec_fd, Nerror_fd, '.', 'markersize', 15, 'DisplayName', 'Finite Difference Error');
semilogy(Nvec_fd, leading_taylor_error, '--', 'DisplayName', '$\frac{1}{30} h^4 \|u^{(5)}\|$');

xlabel('N');
ylabel('Error');
title('Convergence of finite differences');
legend('Location', 'southwest', 'interpreter', 'latex');

% Plot difference between Taylor series and finite difference
%
subplot(2, 1, 2);

hold on;
grid on;
loglog(Nvec_fd, abs(Nerror_fd - leading_taylor_error), '.', 'markersize', 15, 'DisplayName', 'Taylor Remainder [$\mathcal{O} (h^6)$]');
semilogy(Nvec_fd, Nvec_fd.^(-5), '--', 'DisplayName', '$N^{-5}$');

xlabel('N');
ylabel('Error');
title('Taylor Series Remainder');
legend('Location', 'southwest', 'interpreter', 'latex');

print('-dpng', 'problem_3')