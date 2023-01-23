u = @(x) exp(sin(x));
uprime = @(x) cos(x) .* u(x);

fd_diff_matrix(u, uprime, 'problem_1_fd')