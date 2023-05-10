% Iterated_MLSApproxApprox2D
% Script that performs iterated approximate MLS approximation
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2);
  D = 64/9;  % Parameter for basis function
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  neval = 40;
  N = 289; gridtype = 'h';
  % Convert D to epsilon for use with basis function definition
  h = 1/(sqrt(N)-1); ep = 1/(sqrt(D)*h);
  % Number of levels for multilevel iteration
  maxlevel = 10000;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  respoints = dsites; ctrs = dsites;
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute exact solution
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute evaluation matrix directly based on the distances
  % between the evaluation points and centers
  DM = DistanceMatrix(epoints,ctrs);
  EM = rbf(ep,DM)/(pi*D);
  % Compute - for all levels - evaluation matrices for
  % residuals directly based on the distances between the
  % next finer points (respoints) and centers
  DM = DistanceMatrix(respoints,ctrs);
  RM = rbf(ep,DM)/(pi*D);
  Pf = zeros(neval^2,1);  %  initialize
  % Create vector of function values (initial residual),
  rhs = testfunction(dsites(:,1),dsites(:,2));
  for level=1:maxlevel
     % Update on evaluation points
     % (for error computation and plotting)
     Pf = Pf + EM*rhs;
     % Compute new residual on data sites
     rhs = rhs - RM*rhs;
     % Compute errors on evaluation grid
     maxerr(level) = norm(Pf-exact,inf);
     rms_err(level) = norm(Pf-exact)/neval;
  end
  figure; semilogy(1:maxlevel,maxerr,'b',1:maxlevel,rms_err,'r');

