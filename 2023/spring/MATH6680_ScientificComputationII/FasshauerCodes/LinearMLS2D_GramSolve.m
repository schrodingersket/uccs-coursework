% LinearMLS2D_GramSolve
% Script that performs MLS approximation with linear reproduction
% by solving local Gram systems
% Calls on: DistanceMatrix, PlotSurf, PlotError2D
  rbf = @(e,r) exp(-(e*r).^2); ep = 5.5;
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  N = 1089; gridtype = 'u';
  neval = 40;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;
  % Create vector of function (data) values
  f = testfunction(dsites(:,1),dsites(:,2));
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute MLS approximation with shifted basis polynomials
  Pf = zeros(neval^2,1);
  for j=1:neval^2
     P = [ones(N,1) dsites-repmat(epoints(j,:),N,1)];
     W = diag(rbf(ep,DistanceMatrix(epoints(j,:),ctrs)));
     c = (P'*W*P)\(P'*W*f);
     Pf(j) = c(1);
  end
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute errors on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  % Plot interpolant
  caption = 'MLS approximation false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,[160,20],caption);
  % Plot absolute error
  caption = 'Error for MLS approximation.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,[160,20],caption)
