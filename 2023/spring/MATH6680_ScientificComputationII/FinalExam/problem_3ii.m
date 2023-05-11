clf;
clear all;


% Gaussian RBF
%
rbf = @(e, r) exp(-(e*r).^2);
d2rbf = @(e, r) 2 * e^2 * (2 * (e*r).^2 - 1) .* rbf(e, r);

N = 24;

% Generate RBF differentiation matrix for second derivative
%
[D2,x] = D2RBF(N,rbf,d2rbf); 
y = x;

[xx,yy] = meshgrid(x,y);

xx = xx(:); 
yy = yy(:);

I = eye(N+1);
L = kron(I,D2) + kron(D2,I);

% Impose boundary conditions by replacing appropriate rows of L
%
b = find(abs(xx)==1 | abs(yy)==1);            % boundary pts
L(b,:) = zeros(4*N,(N+1)^2); 
L(b,b) = eye(4*N);

% Forcing function for Poisson equation
%
f = -exp(-(yy.^2 + xx.^2));
f(b) = zeros(4*N,1);

% Solve for u and reshape single vector to 2D grid
%
u = L\f;
uu = reshape(u,N+1,N+1);

% Plot solution
%
[xx,yy] = meshgrid(x,y);
[xxx,yyy] = meshgrid(-1:.0333:1,-1:.0333:1);
uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');

surf(xxx,yyy,uuu)
xlabel('x')
ylabel('y')
zlabel('u')
