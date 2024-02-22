% LOOCV2Dmin
% Script that performs leave-one-out cross-validation
% (Rippa's method) to find a good epsilon for 2D RBF interpolation
% with the help of Matlab's fminbnd
% Calls on: DistanceMatrix
% Requires: CostEpsilon
  rbf = @(e,r) exp(-(e*r).^2);   % Gaussian RBF
  % Parameters for shape parameter optimization below
  mine = 0; maxe = 20;
  % Number and type of data points
  N = 81; gridtype = 'u';
  % Define test function
  testfunction = @(x,y) sinc(x).*sinc(y);
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;   % centers coincide with data sites
  % Create right-hand side vector, i.e.,
  % evaluate the test function at the data points.
  rhs = testfunction(dsites(:,1),dsites(:,2));
  % Compute distance matrix between the data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  [ep,fval] = fminbnd(@(ep) CostEpsilon(ep,DM_data,rbf,rhs),...
                      mine,maxe);
  fprintf('Smallest maximum norm: %e\n', fval)
  fprintf('at epsilon = %f\n', ep)
 