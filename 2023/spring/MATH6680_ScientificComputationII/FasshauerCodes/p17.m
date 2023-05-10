% p17
% Script that solves Helmholtz equation
% u_xx + u_yy + (k^2)u = f      on [-1,1]x[-1,1]
% We replace the Chebyshev method by an RBF-PS method
% and explicitly enforce the boundary conditions
% Calls on: D2RBF
  % Gaussian RBF basic function
  rbf = @(e,r) exp(-(e*r).^2);
  d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
  N = 24;
  [D2,x] = D2RBF(N,rbf,d2rbf); y = x;
  [xx,yy] = meshgrid(x,y);
  xx = xx(:); yy = yy(:);
  I = eye(N+1);
  k = 9;
  L = kron(I,D2) + kron(D2,I) + k^2*eye((N+1)^2);
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
