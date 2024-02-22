% RBFInterpolation2D
% Script that performs basic 2D RBF interpolation
% Calls on: DistanceMatrix
  % Define the Gaussian RBF and shape parameter
  rbf = @(e,r) exp(-(e*r).^2); ep = 21.1;
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  N = 1089; gridtype = 'h';
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;
  neval = 40; grid = linspace(0,1,neval);
  [xe,ye] = meshgrid(grid); epoints = [xe(:) ye(:)];
  % Evaluate the test function at the data points
  rhs = testfunction(dsites(:,1),dsites(:,2));
  % Compute distance matrix between the data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  % Compute interpolation matrix
  IM = rbf(ep,DM_data);
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  % Compute evaluation matrix
  EM = rbf(ep,DM_eval);
  % Compute RBF interpolant
  % (evaluation matrix * solution of interpolation system)
  Pf = EM * (IM\rhs);
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute errors on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  % Plot interpolant
  fview = [160,20]; % for Franke's function
  caption = ['RBF interpolant ',...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  % Plot maximum error
  caption = 'Maximum error for RBF interpolant.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)

