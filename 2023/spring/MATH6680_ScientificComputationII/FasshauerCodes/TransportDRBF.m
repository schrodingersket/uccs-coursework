% TransportDRBF
% Script that solves constant coefficient wave equation
% u_t + c*u_x = 0, using RBF-PS approach
% Calls on: DRBF
  rbf = @(e,r) exp(-(e*r).^2);     % Gaussian RBF
  dxrbf = @(e,r,dx) -2*dx*e^2.*exp(-(e*r).^2);
  N = 20;
  [D,x] = DRBF(N,rbf,dxrbf);
  x = flipud(x); dt = 0.001; t = 0; c = -1;
  v = 64*(-x).^3.*(1+x).^3;
  v(find(x>0)) = zeros(length(find(x>0)),1);
  % Time-stepping by explicit Euler formula:
  tmax = 1; tplot = .02; plotgap = round(tplot/dt);
  dt = tplot/plotgap; nplots = round(tmax/tplot);
  data = [v'; zeros(nplots,N+1)]; tdata = t;
  I = eye(size(D));
  for i = 1:nplots
     for n = 1:plotgap
        t = t+dt;
        vv = v(end-1);
        v = v - dt*c*(D*v);     % explicit Euler
        v(1) = 0; v(end) = vv;
     end
     data(i+1,:) = v'; tdata = [tdata; t];
  end
  surf(x,tdata,data), view(10,70), colormap('default');
  axis([-1 1 0 tmax 0 1]), ylabel t, zlabel u, grid off
  % exact solution and error
  xx = linspace(-1,1,101);
  vone = 64*(1-xx).^3.*xx.^3;
  vone(find(xx<0)) = zeros(length(find(xx<0)),1);
  w = interp1(x,v,xx);
  maxErr = norm(w-vone,inf)
