% p17_2D
% Script that solves Helmholtz equation
% u_xx + u_yy + (k^2)u = f   on [-1,1]x[-1,1]
% We replace the Chebyshev method by an RBF-PS method,
% explicitly enforce the boundary conditions, and
% use a 2-D implementation of the Laplacian
% Calls on: LRBF
  % Wendland C6 RBF basic function
  rbf = @(e,r) max(1-e*r,0).^8.*(32*(e*r).^3+25*(e*r).^2+8*e*r+1);
  Lrbf = @(e,r) 44*e^2*max(1-e*r,0).^6.*...
                (88*(e*r).^3+3*(e*r).^2-6*e*r-1);
  [L,x,y] = LRBF(N,rbf,Lrbf);
  [xx,yy] = meshgrid(x,y);
  xx = xx(:); yy = yy(:);
  k = 9;
  L = L + k^2*eye((N+1)^2);
  % Impose boundary conditions by replacing appropriate rows of L
  b = find(abs(xx)==1 | abs(yy)==1);            % boundary pts
  L(b,:) = zeros(4*N,(N+1)^2); L(b,b) = eye(4*N);
  f = exp(-10*((yy-1).^2+(xx-.5).^2));
  f(b) = zeros(4*N,1);
  % Solve for u, reshape to 2D grid, and plot:
  u = L\f;
  uu = reshape(u,N+1,N+1);
  [xx,yy] = meshgrid(x,y);
  [xxx,yyy] = meshgrid(-1:.0333:1,-1:.0333:1);
  uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
  figure, clf, surf(xxx,yyy,uuu),
  xlabel x, ylabel y, zlabel u
  text(.2,1,.022,sprintf('u(0,0) = %13.11f',uu(N/2+1,N/2+1)))
