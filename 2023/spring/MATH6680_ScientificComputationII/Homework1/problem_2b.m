u = @(x) exp(sin(x) .* abs(sin(x)));
uprime = @(x) (sin(x) .* sin(2*x) .* u(x)) ./ abs(sin(x));

fd_diff_matrix(u, uprime, 'problem_2b_fd')
spectral_diff_matrix(u, uprime, 'problem_2b_spectral')