% HermiteLaplaceMixedBCTref_2D
% Script that performs Hermite collocation for 2D Laplace equation
% Note: Prog 36 in Trefethen (2000), exact solution not provided
% Calls on: DistanceMatrix
  % IMQ RBF and its Laplacian
  rbf = @(e,r) 1./sqrt(1+(e*r).^2); ep = 3;
  Lrbf = @(e,r) e^2*((e*r).^2-2)./(1+(e*r).^2).^(5/2);
  L2rbf = @(e,r) 3*e^4*(3*(e*r).^4-24*(e*r).^2+8)./...
                       (1+(e*r).^2).^(9/2);
  %  Laplacian for test problem
  Lu = @(x,y) zeros(size(x));
  % Number and type of collocation points
  N = 289; gridtype = 'u';
  neval = 41;
  % Load (interior) collocation points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name);
  intdata = 2*dsites-1;
  % Additional boundary collocation points
  sn = sqrt(N); bdylin = linspace(-1,1,sn)';
  bdy1 = ones(sn-1,1);
  bdydata = [bdylin(1:end-1) -bdy1; bdy1 bdylin(1:end-1); ...
       flipud(bdylin(2:end)) bdy1; -bdy1 flipud(bdylin(2:end))];
  % Create additional boundary centers OUTSIDE the domain
  h = 2/(sn-1); bdylin = (-1+h:h:1-h)';
  bdy0 = repmat(-1-h,sn-2,1); bdy1 = repmat(1+h,sn-2,1);
  bdyctrs = [-1-h -1-h; bdylin bdy0; 1+h -1-h; bdy1 bdylin; ...
       1+h 1+h; flipud(bdylin) bdy1; -1-h 1+h; bdy0 flipud(bdylin)];
  intctrs = intdata;
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(-1,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute evaluation matrix
  DM_inteval = DistanceMatrix(epoints,intctrs);
  LEM = Lrbf(ep,DM_inteval);
  DM_bdyeval = DistanceMatrix(epoints,bdyctrs);
  BEM = rbf(ep,DM_bdyeval);
  EM = [LEM BEM];
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
  rhs = [Lu(intdata(:,1),intdata(:,2)); zeros(sn-1,1); ...
         0.2*sin(3*pi*bdydata(sn:2*sn-2,2)); zeros((sn-1)/2,1);...
         sin(pi*bdydata((5*sn-3)/2:3*sn-3,1)).^4; zeros(sn-1,1)];
  % Compute RBF solution
  Pf = EM * (CM\rhs);
  surf(xe,ye,reshape(Pf,neval,neval));
  view(-20,45),  axis([-1 1 -1 1 -.2 1]);
  text(0,.8,.5,sprintf('u(0,0) = %12.10f',Pf(841)))
