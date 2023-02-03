function [Nvec, Nerror, T] = fd_diff_matrix(ufn, uprimefn)
  % p1.m - convergence of fourth-order finite differences
  % For various N, set up grid in [-pi,pi] and function u(x):
  Nvec = 2.^(3:16);

  Nerror = zeros(size(Nvec));
  T = zeros(size(Nvec));

  for i = 1:length(Nvec)
    N = Nvec(i);
    h = 2*pi/N;
    x = -pi + (1:N)' * h;
    u = ufn(x);
    uprime = uprimefn(x);

    % Construct sparse 4th-order differentiation matrix:
    %
    tStart = tic;
    e = ones(N, 1);
    D = sparse(1:N, [2:N 1], 2*e/3, N, N) - sparse(1:N, [3:N 1 2], e/12, N, N);
    D = (D - D')/h;

    % Store finite difference computation time and error
    %
    Nerror(i) = norm(D * u - uprime, inf);
    T(i) = toc(tStart);
  end
end
