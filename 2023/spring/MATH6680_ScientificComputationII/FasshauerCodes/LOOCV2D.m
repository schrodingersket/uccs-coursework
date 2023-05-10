% LOOCV2D
% Script that performs leave-one-out cross-validation
% (Rippa's method) to find a good epsilon for 2D RBF interpolation
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2);   % Gaussian RBF
  % Parameters for shape parameter loop below
  mine = 0; maxe = 20; ne = 500;
  ep = linspace(mine,maxe,ne);
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
  for i=1:length(ep)
     % Compute interpolation matrix
     IM = rbf(ep(i),DM_data);
     % Compute error function (i.e., "cost" of epsilon)
     invIM = pinv(IM);
     EF = (invIM*rhs)./diag(invIM);
     % Compute maximum norm of EF
     maxEF(i) = norm(EF(:),inf);
  end
  fprintf('Smallest maximum norm: %e\n', min(maxEF))
  fprintf('at epsilon = %f\n',ep(maxEF==min(maxEF)))
  % Plot cost function norm
  figure; semilogy(ep,maxEF,'b');
