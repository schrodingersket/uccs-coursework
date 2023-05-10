% RBFInterpolation2DtpsH
% Script that performs 2D TPS interpolation with homogeneous kernel
% Calls on: tpsH
  function RBFInterpolation2DtpsH
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  N = 25; gridtype = 'u';
  a = 1e9;
  ppoints = a*[0 0; 1 0; 0 1];
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  % Remove (0,0), (1,0), (0,1) to work with C matrix
  remove = [find(dsites(:,1)==0 & dsites(:,2)==0);...
       find(dsites(:,1)==1 & dsites(:,2)==0);...
       find(dsites(:,1)==0 & dsites(:,2)==1)];
  dsites(remove,:) = [];
  % Scale problem to square [0,a]^2
  dsites = a*dsites;
  % Let centers coincide with data sites
  ctrs=dsites;
  neval = 40; grid = linspace(0,a,neval);
  [xe,ye] = meshgrid(grid); epoints = [xe(:) ye(:)];
  % Create right-hand side for homogeneous problem
  DP = [p1(dsites/a) p2(dsites/a) p3(dsites/a)]';
  d = testfunction(ppoints(:,1)/a,ppoints(:,2)/a);
  rhs = testfunction(dsites(:,1)/a,dsites(:,2)/a) - DP'*d;
  % Compute interpolation matrix for the special case of TPS
  % native space kernel (no need to add polynomials)
  IM = tpsH(dsites,ctrs,a);
  % Compute condition number of interpolation matrix
  fprintf('l2-condition         : %e\n', cond(IM))
  % Compute evaluation matrix
  EM = tpsH(epoints,ctrs,a);
  EP = [p1(epoints/a) p2(epoints/a) p3(epoints/a)];
  EM = [EM EP];
  % Compute RBF interpolant
  Pf = EM * [(IM\rhs); d];
  % Compute exact solution
  exact = testfunction(epoints(:,1)/a,epoints(:,2)/a);
  % Compute errors on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  fview = [160,20];  % viewing angles for plot
  caption = ['Homogeneous TPS interpolant '...
      'false colored by maximum error.'];
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for homogeneous TPS interpolant.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
  return
  % The cardinal polynomials
  function w = p1(z)
  w = 1 - z(:,1) - z(:,2);
  return
  function w = p2(z)
  w = z(:,1);
  return
  function w = p3(z)
  w = z(:,2);
  return
