% RBFKnotInsert2D
% Script that performs 2D RBF least squares approximation
% via knot insertion
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2); ep = 5.5;
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  N = 289; gridtype = 'h';
  M = 1;  % Number of initial centers
  neval = 40;
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  tol = 1e-5;  % Tolerance; stopping criterion
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  % Take first M "data sites" as centers
  ctrs = dsites(1:M,:);
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Create right-hand side vector, i.e.,
  % evaluate the test function at the data points.
  rhs = testfunction(dsites(:,1),dsites(:,2));
  rms_res = 999999;
  while (rms_res > tol)
     % Compute least squares fit
     DM_data = DistanceMatrix(dsites,ctrs);
     CM = rbf(ep,DM_data);
     coef = CM\rhs;
     % Compute residual
     residual = abs(CM*coef - rhs);
     [sresidual,idx] = sort(residual);
     lres = length(residual);
     rms_res = norm(residual)/sqrt(lres);
     % Add point(s)
     if (rms_res > tol)
        addpoint = idx(lres);   % This is the point we add
        % If already used, try next point
        while any(ismember(ctrs,dsites(addpoint,:),'rows'))
           lres = lres-1;  addpoint = idx(lres);
        end
        ctrs = [ctrs; dsites(addpoint,:)];
     end
  end
  % Compute evaluation matrix
  DM_eval = DistanceMatrix(epoints,ctrs);
  EM = rbf(ep,DM_eval);
  Pf = EM*coef;   % Compute RBF least squares approximation
  % Compute maximum error on evaluation grid
  maxerr = max(abs(Pf - exact));  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  figure; % Plot data sites and centers
  plot(dsites(:,1),dsites(:,2),'bo',ctrs(:,1),ctrs(:,2),'r+');
  caption = sprintf('%d data sites and %d centers',N,size(ctrs,1));
  title(caption);
  caption = 'RBF approximant false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,[160,20],caption);
