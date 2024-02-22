% Shepard_CS
% Script that performs Shepard approximation for arbitrary
% space dimensions s using sparse matrices
% Calls on: DistanceMatrixCSRBF, MakeSDGrid, testfunction
% Uses:     haltonseq (written by Daniel Dougherty from
%           Matlab Central File Exchange)
  % Wendland C2 weight function
  rbf = @(e,r) r.^4.*(5*spones(r)-4*r);
  s = 6;   % Space dimension s
  % Number of Halton data points
  k = 3; N = (2^k+1)^s;
  ep = 2^k+1; % Scale parameter for basis function
  neval = 4; M = neval^s;
  % Compute data sites as Halton points
  dsites = haltonseq(N,s);
  ctrs = dsites;
  % Create vector of function (data) values
  f = testfunction(s,dsites);
  % Create neval^s equally spaced evaluation locations in
  % the s-dimensional unit cube
  epoints = MakeSDGrid(s,neval);
  % Compute evaluation matrix, i.e.,
  % matrix of values of generating functions
  DM_eval = DistanceMatrixCSRBF(epoints,ctrs,ep);
  EM = rbf(ep,DM_eval);
  % Shepard scaling
  EM = spdiags(1./(EM*ones(N,1)),0,M,M)*EM;
  % Compute quasi-interpolant
  Pf = EM*f;
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunction(s,epoints);
  % Compute errors on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/sqrt(M);
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
