% RBFApproximation2D
% Script that performs basic 2D RBF least squares approximation
% Calls on: DistanceMatrix, PlotSurf, PlotError2D
  rbf = @(e,r) exp(-(e*r).^2); ep = 1;
  testfunction = @(x,y) sinc(x).*sinc(y);
  N = 1089; gridtype = 'h';
  M = 81; grid2type = 'u';
  neval = 40;
  % Load centers
  name = sprintf('Data2D_%d%s',M,grid2type); load(name)
  ctrs = dsites;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  % Compute distance matrix between data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  % Build collocation matrix
  CM = rbf(ep,DM_data);
  % Create right-hand side vector, i.e.,
  % evaluate the test function at the data points.
  rhs = testfunction(dsites(:,1),dsites(:,2));
  % Create neval-by-neval equally spaced evaluation
  % locations in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  EM = rbf(ep,DM_eval);
  % Compute RBF least squares approximation
  Pf = EM * (CM\rhs);
  % Compute exact solution, i.e., evaluate test
  % function on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute maximum error on evaluation grid
  maxerr = norm(Pf-exact,inf);
  % Plots
  figure; fview = [100,30]; % viewing angles for plot
  caption = sprintf('%d data sites and %d centers',N,M);
  title(caption);
  plot(dsites(:,1),dsites(:,2),'bo',ctrs(:,1),ctrs(:,2),'r+');
  caption = 'RBF approximant false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for RBF approximant.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
