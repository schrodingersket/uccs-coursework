% KansaEllipticVC_2D
% Script that performs Kansa collocation for 2D elliptic PDE
% with variable coefficients
% Calls on: DistanceMatrix, DifferenceMatrix
  % IMQ RBF and its derivatives
  rbf = @(e,r) 1./sqrt(1+(e*r).^2); ep = 3;
  dxrbf = @(e,r,dx) -dx*e^2./(1+(e*r).^2).^(3/2);
  dyrbf = @(e,r,dy) -dy*e^2./(1+(e*r).^2).^(3/2);
  dxxrbf = @(e,r,dx) e^2*(3*(e*dx).^2-1-(e*r).^2)./...
                         (1+(e*r).^2).^(5/2);
  dyyrbf = @(e,r,dy) e^2*(3*(e*dy).^2-1-(e*r).^2)./...
                         (1+(e*r).^2).^(5/2);
  % Test problem input (right-hand side, coefficients)
  u = @(x,y) 16*x.*(1-x).*y.*(1-y);
  Lu = @(x,y) -16*x.*exp(x-y).*(1-x).*(3-2*y)+...
             32*y.*(1-y).*(3*x.^2+y.^2-x-2);
  a = @(x,y) 2-x.^2-y.^2; ax = @(x,y) -2*x;
  b = @(x,y) exp(x-y); by = @(x,y)-exp(x-y);
  N = 289; gridtype = 'h';
  neval = 40;
  % Load (interior) collocation points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name);
  intdata = dsites;
  % Additional boundary collocation points
  sn = sqrt(N); bdylin = linspace(0,1,sn)';
  bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
  bdydata = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1);...
       flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
  % Create additional boundary centers OUTSIDE the domain
  h = 1/(sn-1); bdylin = (h:h:1-h)';
  bdy0 = -h*ones(sn-2,1); bdy1 = (1+h)*ones(sn-2,1);
  bdyctrs = [-h -h; bdylin bdy0; 1+h -h; bdy1 bdylin;...
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
  dx_intdata = Differencematrix(intdata(:,1),ctrs(:,1));
  dy_intdata = Differencematrix(intdata(:,2),ctrs(:,2));
  LCM = diag(ax(intdata(:,1))) * ...
        dxrbf(ep,DM_intdata,dx_intdata) + ...
        diag(a(intdata(:,1),intdata(:,2))) * ...
        dxxrbf(ep,DM_intdata,dx_intdata) + ...
        diag(by(intdata(:,1),intdata(:,2))) * ...
        dyrbf(ep,DM_intdata,dy_intdata) + ...
        diag(b(intdata(:,1),intdata(:,2))) * ...
        dyyrbf(ep,DM_intdata,dy_intdata);
  BCM = rbf(ep,DM_bdydata);
  CM = [LCM; BCM];
  % Create right-hand side
  rhs = [Lu(intdata(:,1),intdata(:,2)); zeros(4*(sn-1),1)];
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
