% TPS_RidgeRegression2D
% Script that performs 2D TPS-RBF approximation with reproduction of
% linear functions and smoothing via ridge regression
% Calls on: tps, DistanceMatrix
  % Use TPS (defined in tps.m, see Appendix C)
  rbf = @tps; ep = 1;
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  omega = 1;   % Smoothing parameter
  N = 1089; gridtype = 'h';
  neval = 40;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;
  % Compute distance matrix between data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  % Create noisy right-hand side vector
  rhs = testfunction(dsites(:,1),dsites(:,2));
  randn('state',3);
  rhs = rhs + 0.03*randn(size(rhs));
  % Add zeros for 2D linear reproduction
  rhs = [rhs; zeros(3,1)];
  % Compute interpolation matrix and add diagonal regularization
  IM = rbf(ep,DM_data);
  IM = IM + eye(size(IM))/(2*omega);
  % Add extra columns and rows for linear reproduction
  PM = [ones(N,1) dsites]; IM = [IM PM; [PM' zeros(3,3)]];
  fprintf('Condition number estimate: %e\n',condest(IM))
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  % Compute evaluation matrix and add columns for linear precision
  EM = rbf(ep,DM_eval);
  PM = [ones(neval^2,1) epoints]; EM = [EM PM];
  % Compute RBF interpolant
  Pf = EM * (IM\rhs);
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute maximum error on evaluation grid
  maxerr = norm(Pf-exact,inf);
  maxerr = max(abs(Pf - exact));
  rms_err = norm(Pf-exact)/neval;
  disp(sprintf('RMS error:     %e', rms_err))
  disp(sprintf('Maximum error: %e', maxerr))
  % Plots
  fview = [160,20]; % viewing angles for plot
  caption = ['TPS ridge regression fit ',...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for TPS rigde regression fit.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
