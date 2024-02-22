% p9.m - polynomial interpolation in equispaced and Chebyshev pts
NN = 1:12;
xx = -1.01:.005:1.01; clf

equispace_idx = 1;
chebyshev_idx = 2;

minimax_error = @(N) 1./(factorial(N + 1) .* 2.^N);
taylor_error = @(N) ((N + 2)./(N + 1)) ./ factorial(N + 1);

polyfit_err = zeros(2, length(NN));

for n = 1:length(NN)
  N = NN(n);
  for i = 1:2
    if i==equispace_idx; x = -1 + 2*(0:N)/N; end
    if i==chebyshev_idx; x = cos(pi*(0:N)/N); end
    u = exp(x);
    uu = exp(xx);
    p = polyfit(x,u,N); % interpolation
    pp = polyval(p,xx); % evaluation of interpolant
    
    polyfit_err(i, n) = norm(uu-pp,inf);
  end
end

hold on;

semilogy(NN, taylor_error(NN), 'k--', 'linewidth', 1, 'DisplayName', 'Truncated Taylor Series');
semilogy(NN, polyfit_err(1, :), 'k', 'linewidth', 1, 'DisplayName', 'Equispaced Interpolation');

semilogy(NN, minimax_error(NN), 'b--', 'linewidth', 1, 'DisplayName', 'Minimax Error');
semilogy(NN, polyfit_err(2, :), 'b', 'linewidth', 1, 'DisplayName', 'Chebyshev Interpolation');

legend;
title('N vs. Error');

xlabel('N');
ylabel('Error');

print('problem_1b', '-dpng');
