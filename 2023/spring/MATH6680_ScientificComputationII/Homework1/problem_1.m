u = @(x) exp(sin(x));
uprime = @(x) cos(x) .* u(x);

[Nvec, Nerror, T] = fd_diff_matrix(u, uprime);

clf;
subplot('position', [.1 .4 .8 .5]);

subplot(2, 1, 1);
hold on;
grid on;

loglog(Nvec, Nerror, '.', 'markersize', 15, 'DisplayName', 'Finite Difference Error')
semilogy(Nvec, Nvec.^(-4), '--', 'DisplayName', 'N^{-4}')

xlabel('N')
ylabel('Error')
title('Convergence of finite differences')
legend('Location', 'southwest');

% Plot partitions vs. computation time
%
subplot(2, 1, 2)

plot(Nvec, T)

xlabel('N')
ylabel('Time')

print('-dpng', 'problem_1_fd')