% KansaLaplaceMixedBC_2D
% Script that performs Kansa collocation for 2D Laplace equation
% with mixed BCs
% Calls on: DistanceMatrix, DifferenceMatrix
  % IMQ RBF and its Laplacian
  rbf = @(e,r) 1./sqrt(1+(e*r).^2); ep = 3;
  dyrbf = @(e,r,dy) -dy*e^2./(1+(e*r).^2).^(3/2);
  Lrbf = @(e,r) e^2*((e*r).^2-2)./(1+(e*r).^2).^(5/2);
  % Exact solution and its Laplacian for test problem
  u = @(x,y) 1-0.9*x.^3+0*y;
  Lu = @(x,y) -5.4*x+0*y;
  % Number and type of collocation points
  N = 289; gridtype = 'h';
  neval = 40;
  % Load (interior) collocation points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name);
  intdata = dsites;
  % Additional boundary collocation points
  sn = sqrt(N); bdylin = linspace(0,1,sn)';
  bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
  bdydata = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1); ...
       flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
  % Create additional boundary centers OUTSIDE the domain
  h = 1/(sn-1); bdylin = (h:h:1-h)';
  bdy0 = -h*ones(sn-2,1); bdy1 = (1+h)*ones(sn-2,1);
  bdyctrs = [-h -h; bdylin bdy0; 1+h -h; bdy1 bdylin; ...
      1+h 1+h; flipud(bdylin) bdy1; -h 1+h; bdy0 flipud(bdylin)];
  ctrs = [intdata; bdyctrs];
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute evaluation matrix
  DM_eval = DistanceMatrix(epoints,ctrs);
  EM = rbf(ep,DM_eval);
  exact = u(epoints(:,1),epoints(:,2));
  % Compute blocks for collocation matrix
  DM_intdata = DistanceMatrix(intdata,ctrs);
  DM_bdydata = DistanceMatrix(bdydata,ctrs);
  dy_bdydata = Differencematrix(bdydata(:,2),ctrs(:,2));
  LCM = Lrbf(ep,DM_intdata);
  BCM1 = -dyrbf(ep,DM_bdydata(1:sn-1,:),dy_bdydata(1:sn-1,:));
  BCM2 = rbf(ep,DM_bdydata(sn:2*sn-2,:));
  BCM3 = dyrbf(ep,DM_bdydata(2*sn-1:3*sn-3,:),...
               dy_bdydata(2*sn-1:3*sn-3,:));
  BCM4 = rbf(ep,DM_bdydata(3*sn-2:end,:));
  CM = [LCM; BCM1; BCM2; BCM3; BCM4];
  % Create right-hand side
  rhs = [Lu(intdata(:,1),intdata(:,2)); zeros(sn-1,1); ...
         0.1*ones(sn-1,1); zeros(sn-1,1); ones(sn-1,1)];
  % RBF solution
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
  caption = ['Nonsymmetric RBF solution '...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for nonsymmetric RBF solution.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
