% RBFApproximation2Dlinear
% Script that performs 2D RBF least squares approximation with
% linear reproduction for noisy data
% Calls on: tps, DistanceMatrix
  rbf = @tps; ep = 1;   % defined in tps.m (see Appendix C)
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
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
  CM = rbf(ep,DM_data);   % Collocation matrix
  % Add extra columns and rows for linear reproduction
  PM = [ones(N,1) dsites]; PtM = [ones(M,1) ctrs]';
  CM = [CM PM; [PtM zeros(3,3)]];
  % Create right-hand side vector and add noise
  rhs = testfunction(dsites(:,1),dsites(:,2));
  rhs = rhs + 0.03*randn(size(rhs));
  % Add zeros for linear (2D) reproduction
  rhs = [rhs; zeros(3,1)];
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  EM = rbf(ep,DM_eval);   % Evaluation matrix
  % Add columns for linear reproduction
  PM = [ones(neval^2,1) epoints]; EM = [EM PM];
  % Compute RBF least squares approximation
  Pf = EM * (CM\rhs);
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute maximum error on evaluation grid
  maxerr = norm(Pf-exact,inf);
  % Plots
  figure; fview = [160,20]; % viewing angles for plot
  caption = sprintf('%d data sites and %d centers',N,M);
  title(caption);
  plot(dsites(:,1),dsites(:,2),'bo',ctrs(:,1),ctrs(:,2),'r+');
  caption = 'RBF approximant false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for RBF approximant.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
