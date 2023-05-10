% RBFGalerkin2D
% Script that performs Galerkin solution of 2D Helmholtz equation
% - u_xx - u_yy + u = f
% Calls on: DistanceMatrix, PlotSurf, PlotError2D
  % Definition of the RBF and its gradient, Wendland C2
  rbf = @(e,x,y,xi,yi) max(1-e*sqrt((x-xi).^2+(y-yi).^2),0).^4.*...
                     (4*e*sqrt((x-xi).^2+(y-yi).^2)+1);
  dxrbf = @(e,x,y,xi,yi) -20*(x-xi)*e^2.*...
                      max(1-e*sqrt((x-xi).^2+(y-yi).^2),0).^3;
  dyrbf = @(e,x,y,xi,yi) -20*(y-yi)*e^2.*...
                      max(1-e*sqrt((x-xi).^2+(y-yi).^2),0).^3;
  evalrbf = @(e,r) max(1-e*r,0).^4.*(4*e*r+1);
  % Products for integration
  rp = @(e,x,y,xi,yi,xj,yj) rbf(e,x,y,xi,yi).*rbf(e,x,y,xj,yj);
  gp = @(e,x,y,xi,yi,xj,yj) dxrbf(e,x,y,xi,yi).*...
         dxrbf(e,x,y,xj,yj)+dyrbf(e,x,y,xi,yi).*dyrbf(e,x,y,xj,yj);
  % Parameter for basis function
  ep = .7;
  % Right-hand side function for Helmholtz equation
  f = @(x,y) cos(pi*x).*cos(pi*y);
  % Exact solution
  u = @(x,y) cos(pi*x).*cos(pi*y)/(2*pi^2+1);
  % Number and type of centers:
  N = 25; gridtype = 'u';
  % Resolution of evaluation grid for errors and plotting
  neval = 40;
  % Load data points
  name = sprintf('Data2D_%d%s', N,gridtype); load(name)
  % Shift centers to the square [-1,1]^2
  ctrs = 2*dsites-1;
  % Build stiffness matrix and right-hand side
  A = zeros(N,N); rhs = zeros(N,1);
  for i=1:N
     for j=1:i
        A(i,j) = dblquad(@(x,y) gp(ep,x,y,ctrs(i,1),ctrs(i,2),...
                            ctrs(j,1),ctrs(j,2)),-1,1,-1,1) + ...
             dblquad(@(x,y) rp(ep,x,y,ctrs(i,1),ctrs(i,2),...
                               ctrs(j,1),ctrs(j,2)),-1,1,-1,1);
     end
     rhs(i) = dblquad(@(x,y) f(x,y).*...
                  rbf(ep,x,y,ctrs(i,1),ctrs(i,2)),-1,1,-1,1);
  end
  % Make matrix symmetric
  A = A + A' - diag(diag(A));
  % Solve linear system, i.e., compute expansion coefficients
  c = A\rhs;
  % Evaluation
  grid = linspace(-1,1,neval); [xe, ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  exact = u(epoints(:,1),epoints(:,2));
  DM_eval = DistanceMatrix(epoints,ctrs);
  EM = evalrbf(ep,DM_eval);
  Pf = EM * c;
  % Compute maximum error on evaluation grid
  maxerr = norm(Pf-exact,inf); rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  % Plot approximate solution
  fview = [-30,30]; % viewing angle for plot
  caption = 'Approximate solution';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
