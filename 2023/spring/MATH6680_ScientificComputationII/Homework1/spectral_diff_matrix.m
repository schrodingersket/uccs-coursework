function [Nvec, Nerror, T] = spectral_diff_matrix(ufn, uprimefn, iterations)
  % p1.m - convergence of fourth-order finite differences
  % For various N, set up grid in [-pi,pi] and function u(x):
  Nvec = 6:2:100;

  Nerror = zeros(size(Nvec));
  T = zeros(size(Nvec));

  for i = 1:length(Nvec)
    N = Nvec(i);
    h = 2*pi/N;
    x = -pi + (1:N)' * h;
    u = ufn(x);
    uprime = uprimefn(x);

    % Construct spectral differentiation matrix:
    %
    tStart = tic;
    column = [0 .5*(-1).^(1:N-1) .* cot((1:N-1) * h/2)];
    D = toeplitz(column, column([1 N:-1:2]));

    % Plot max(abs(D*u-uprime)):
    %
    Nerror(i) = norm(D * u - uprime, inf);
    T(i) = toc(tStart);
  end

end
