% RBFInterpolation2Dlinear
% Script that performs 2D RBF interpolation with reproduction of
% linear functions
% Calls on: DistanceMatrix
  % Define the Gaussian RBF and shape parameter
  rbf = @(e,r) exp(-(e*r).^2); ep = 6;
  % Define linear test function
  testfunction = @(x,y) (x+y)/2;
  % Number and type of data points
  N = 9; gridtype = 'u';
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;
  neval = 40; M = neval^2; grid = linspace(0,1,neval);
  [xe,ye] = meshgrid(grid); epoints = [xe(:) ye(:)];
  % Evaluate the test function at the data points.
  rhs = testfunction(dsites(:,1),dsites(:,2));
  % Add zeros for linear (2D) reproduction
  rhs = [rhs; zeros(3,1)];
  % Compute distance matrix between the data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  % Compute interpolation matrix
  IM = rbf(ep,DM_data);
  % Define 3-column matrix P for linear reproduction
  PM = [ones(N,1) dsites];
  % Augment interpolation matrix
  IM = [IM PM; [PM' zeros(3,3)]];
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  % Compute evaluation matrix
  EM = rbf(ep,DM_eval);
  % Add column for constant reproduction
  PM = [ones(M,1) epoints]; EM = [EM PM];
  % Compute RBF interpolant
  % (evaluation matrix * solution of interpolation system)
  Pf = EM * (IM\rhs);
  % Compute maximum error on evaluation grid
  exact = testfunction(epoints(:,1),epoints(:,2));
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  % Plot interpolant
  fview = [-30,30];
  caption = ['RBF interpolant (linear reproduction) ',...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  % Plot maximum error
  caption = 'Maximum error for RBF interpolant.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)

