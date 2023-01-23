u = @(x) exp(sin(x).^2);
uprime = @(x) 2 * sin(x) .* cos(x) .* u(x);

fd_diff_matrix(u, uprime, 'problem_2a_fd')
spectral_diff_matrix(u, uprime, 'problem_2a_spectral')