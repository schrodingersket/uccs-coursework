% Powerfunction2D
% Script that finds "optimal" shape parameter by computing the power
% function for the 2D RBF interpolation approach with varying epsilon
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2);   % Define the Gaussian RBF
  % Parameters for shape parameter loop below
  mine = 0; maxe = 20;
  ne = 500; ep = linspace(mine,maxe,ne);
  % Number and type of data points
  N = 81; gridtype = 'u';
  % Resolution of grid for power function norm computation
  neval = 40; M = neval^2;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;   % centers coincide with data sites
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  % Compute distance matrix between the data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  for i=1:length(ep)
     % Compute interpolation matrix
     IM = rbf(ep(i),DM_data);
     % Compute evaluation matrix
     EM = rbf(ep(i),DM_eval);
     % Compute power function at evaluation points
     invIM = pinv(IM); phi0 = rbf(ep(i),0);
     for j=1:M
        powfun(j) = real(sqrt(phi0-(invIM*EM(j,:)')'*EM(j,:)'));
     end
     % Compute max. norm of power function on evaluation grid
     maxPF(i) = max(powfun);
  end
  fprintf('Smallest maximum norm: %e\n', min(maxPF))
  fprintf('at epsilon = %f\n',ep(maxPF==min(maxPF)))
  fprintf('with cond(A) = %e\n', ...
      condest(rbf(ep(find(maxPF==min(maxPF))),DM_data)))
  % Plot power function norm
  figure; semilogy(ep,maxPF,'b');
