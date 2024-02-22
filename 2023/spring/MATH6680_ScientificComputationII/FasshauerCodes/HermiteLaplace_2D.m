% HermiteLaplace_2D
% Script that performs Hermite collocation for 2D Laplace equation
% Calls on: DistanceMatrix
  % IMQ RBF and its Laplacian and double Laplacian
  rbf = @(e,r) 1./sqrt(1+(e*r).^2); ep = 3;
  Lrbf = @(e,r) e^2*((e*r).^2-2)./(1+(e*r).^2).^(5/2);
  L2rbf = @(e,r) 3*e^4*(3*(e*r).^4-24*(e*r).^2+8)./...
                        (1+(e*r).^2).^(9/2);
  % Exact solution and its Laplacian for test problem
  u = @(x,y) sin(pi*x).*cos(pi*y/2);
  Lu = @(x,y) -1.25*pi^2*sin(pi*x).*cos(pi*y/2);
  % Number and type of collocation points
  N = 289; gridtype = 'u';
  neval = 40;
  % Load (interior) collocation points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name);
  intdata = dsites;
  % Additional (equally spaced) boundary collocation points
  sn = sqrt(N); bdylin = linspace(0,1,sn)';
  bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
  bdydata = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1); ...
       flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
  % Create additional boundary centers OUTSIDE the domain
  h = 1/(sn-1); bdylin = (h:h:1-h)';
  bdy0 = -h*ones(sn-2,1); bdy1 = (1+h)*ones(sn-2,1);
  bdyctrs = [-h -h; bdylin bdy0; 1+h -h; bdy1 bdylin; ...
       1+h 1+h; flipud(bdylin) bdy1; -h 1+h; bdy0 flipud(bdylin)];
  intctrs = intdata;
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute evaluation matrix
  DM_inteval = DistanceMatrix(epoints,intctrs);
  LEM = Lrbf(ep,DM_inteval);
  DM_bdyeval = DistanceMatrix(epoints,bdyctrs);
  BEM = rbf(ep,DM_bdyeval);
  EM = [LEM BEM];
  exact = u(epoints(:,1),epoints(:,2));
  % Compute blocks for collocation matrix
  DM_IIdata = DistanceMatrix(intdata,intctrs);
  LLCM = L2rbf(ep,DM_IIdata);
  DM_IBdata = DistanceMatrix(intdata,bdyctrs);
  LBCM = Lrbf(ep,DM_IBdata);
  DM_BIdata = DistanceMatrix(bdydata,intctrs);
  BLCM = Lrbf(ep,DM_BIdata);
  DM_BBdata = DistanceMatrix(bdydata,bdyctrs);
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
  caption = ['Symmetric RBF solution '...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for symmetric RBF solution.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
