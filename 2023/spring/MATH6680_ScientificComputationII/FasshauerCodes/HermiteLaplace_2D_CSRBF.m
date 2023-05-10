% HermiteLaplace_2D_CSRBF
% Script that performs Hermite collocation for 2D Laplace equation
% with sparse matrices
% Calls on: DistanceMatrixCSRBF
  % Wendland C6 RBF, its Laplacian and double Laplacian
  rbf = @(e,r) r.^8.*(66*spones(r)-154*r+121*r.^2-32*r.^3);
  Lrbf = @(e,r) 44*e^2*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
  L2rbf = @(e,r) 1056*e^4*r.^4.*...
                 (105*spones(r)-483*r+679*r.^2-297*r.^3);
  ep = 0.25;
  % Exact solution and its Laplacian for test problem
  u = @(x,y) sin(pi*x).*cos(pi*y/2);
  Lu = @(x,y) -1.25*pi^2*sin(pi*x).*cos(pi*y/2);
  % Number and type of collocation points
  N = 289; gridtype = 'h';
  neval = 40;
  % Load (interior) collocation points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name);
  intdata = dsites;
  % Additional (equally spaced) boundary collocation points
  sn = sqrt(N); bdylin = linspace(0,1,sn)';
  bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
  bdydata = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1); ...
       flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
  % Let centers coincide with ALL data sites
  bdyctrs = bdydata;
  ctrs = [intdata; bdyctrs];
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute evaluation matrix
  DM_inteval = DistanceMatrixCSRBF(epoints,intdata,ep);
  DM_bdyeval = DistanceMatrixCSRBF(epoints,bdyctrs,ep);
  LEM = Lrbf(ep,DM_inteval);
  BEM = rbf(ep,DM_bdyeval);
  EM = [LEM BEM];
  exact = u(epoints(:,1),epoints(:,2));
  % Compute blocks for collocation matrix
  DM_IIdata = DistanceMatrixCSRBF(intdata,intdata,ep);
  DM_IBdata = DistanceMatrixCSRBF(intdata,bdyctrs,ep);
  DM_BIdata = DistanceMatrixCSRBF(bdydata,intdata,ep);
  DM_BBdata = DistanceMatrixCSRBF(bdydata,bdyctrs,ep);
  LLCM = L2rbf(ep,DM_IIdata);
  LBCM = Lrbf(ep,DM_IBdata);
  BLCM = Lrbf(ep,DM_BIdata);
  BBCM = rbf(ep,DM_BBdata);
  CM = [LLCM LBCM; BLCM BBCM];
  % Create right-hand side
  rhs = [Lu(intdata(:,1),intdata(:,2)); ...
         sin(pi*bdydata(1:sn-1,1)); zeros(3*(sn-1),1)];
  % Compute RBF solution
  Pf = EM * (CM\rhs);
  % Compute maximum error on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  % Plot collocation points and centers
  hold on; plot(intdata(:,1),intdata(:,2),'bo');
  plot(bdydata(:,1),bdydata(:,2),'rx');
  plot(bdyctrs(:,1),bdyctrs(:,2),'gx'); hold off
  fview = [-30,30];  % viewing angles for plot
  caption = ['Symmetric CSRBF solution '...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for symmetric CSRBF solution.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
