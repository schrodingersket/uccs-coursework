% RBFGreedyOnePoint2D
% Script that performs greedy one point algorithm for adaptive
% 2D RBF interpolation
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2);     % Gaussian RBF
  ep = 5.5;  % Parameter for basis function
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  % Number and type of data points
  N = 16641; gridtype = 'h';
  neval = 40; grid = linspace(0,1,neval);
  [xe,ye] = meshgrid(grid); epoints = [xe(:) ye(:)];
  % Tolerance; stopping criterion
  tol = 1e-5; kmax = 1000;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  % Initialize residual and fit
  r_old = testfunction(dsites(:,1),dsites(:,2));
  u_old = 0;
  k = 1; maxres(k) = 999999;
  % Use an (arbitrary) initial point
  ykidx = (N+1)/2; yk(k,:) = dsites(ykidx,:);
  while (maxres(k) > tol && k < kmax)
     % Evaluate basis function at yk
     DM_data = DistanceMatrix(yk(k,:),yk(k,:));
     IM = rbf(ep,DM_data);
     beta = r_old(ykidx)/IM;
     % Compute evaluation matrices for residual and fit
     DM_res = DistanceMatrix(dsites,yk(k,:));
     RM = rbf(ep,DM_res);
     DM_eval = DistanceMatrix(epoints,yk(k,:));
     EM = rbf(ep,DM_eval);
     % Update residual and fit
     r = r_old - beta*RM;    u = u_old + beta*EM;
     % Find new point to add
     [sr,idx] = sort(abs(r));
     maxres(k+1) = sr(end);
     ykidx = idx(end);    yk(k+1,:) = dsites(ykidx,:);
     r_old = r;    u_old = u;
     k = k + 1;
  end
  % Compute exact solution
  exact = testfunction(epoints(:,1),epoints(:,2));
  maxerr = norm(u-exact,inf); rms_err = norm(u-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  fview = [160,20]; % viewing angles for plot
  caption = 'Greedy one-point approximation';
  PlotSurf(xe,ye,u,neval,exact,maxerr,fview,caption);
  caption = 'Absolute error';
  PlotError2D(xe,ye,u,exact,maxerr,neval,fview,caption)
  figure; plot(yk(:,1),yk(:,2),'ro')
  figure; semilogy(maxres,'b');
