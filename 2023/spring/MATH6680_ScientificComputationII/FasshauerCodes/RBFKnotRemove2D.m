% RBFKnotRemove2D
% Script that performs 2D RBF least squares approximation
% via knot removal
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2); ep = 5.5;
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  N = 289; gridtype = 'h';
  M = 289;  % Number of initial centers
  neval = 40;
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  tol = 5e-1;  % Tolerance; stopping criterion
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  % Take first M "data sites" as centers
  ctrs = dsites(1:M,:);
  % Compute exact solution, i.e., evaluate test function
  % on evaluation points
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Create right-hand side vector, i.e.,
  % evaluate the test function at the data points.
  rhs = testfunction(dsites(:,1),dsites(:,2));
  minres = 0;
  while (minres < tol)
     % Compute collocation matrix
     DM_data = DistanceMatrix(dsites,ctrs);
     CM = rbf(ep,DM_data);
     % Compute residual
     invCM = pinv(CM);  EF = (invCM*rhs)./diag(invCM);
     residual = abs(EF);
     [sresidual,idx] = sort(residual);  minres = residual(1);
     % Remove point
     if (minres < tol)
        ctrs = [ctrs(1:idx(1)-1,:); ctrs(idx(1)+1:M,:)];
        M = M-1;
     end
  end
  % Evaluate final least squares fit
  DM_data = DistanceMatrix(dsites,ctrs);
  CM = rbf(ep,DM_data);
  DM_eval = DistanceMatrix(epoints,ctrs);
  EM = rbf(ep,DM_eval);
  Pf = EM*(CM\rhs);
  maxerr = max(abs(Pf - exact));  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  figure;   % Plot data sites and centers
  plot(dsites(:,1),dsites(:,2),'bo',ctrs(:,1),ctrs(:,2),'r+');
  caption = sprintf('%d data sites and %d centers', N, M);
  title(caption);
  caption = 'RBF approximant false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,[160,20],caption);
