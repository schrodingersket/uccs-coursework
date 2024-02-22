% p35
% Script that solves Allen-Cahn equation with boundary condition
% imposed explicitly ("method (II)") (from Trefethen (2000))
% We replace the Chebyshev method by an RBF-PS method
% Calls on: D2RBF
  % Matern cubic as RBF basic function
  rbf = @(e,r) exp(-e*r).*(15+15*e*r+6*(e*r).^2+(e*r).^3);
  d2rbf = @(e,r) e^2*((e*r).^3-3*e*r-3).*exp(-e*r);
  N = 20;
  [D2,x] = D2RBF(N,rbf,d2rbf);
  % Here is the rest of Trefethen's code.
  mu = 0.01; dt = min([.01,50*N^(-4)/mu]);
  t = 0; v = .53*x + .47*sin(-1.5*pi*x);
  % Solve PDE by Euler formula and plot results:
  tmax = 100; tplot = 2; nplots = round(tmax/tplot);
  plotgap = round(tplot/dt); dt = tplot/plotgap;
  xx = -1:.025:1; vv = polyval(polyfit(x,v,N),xx);
  plotdata = [vv; zeros(nplots,length(xx))]; tdata = t;
  for i = 1:nplots
     for n = 1:plotgap
        t = t+dt; v = v + dt*(mu*D2*v + v - v.^3); % Euler
        v(1) = 1 + sin(t/5)^2; v(end) = -1; % BC
     end
     vv = polyval(polyfit(x,v,N),xx);
     plotdata(i+1,:) = vv; tdata = [tdata; t];
  end
  surf(xx,tdata,plotdata), grid on
  axis([-1 1 0 tmax -1 2]), view(-40,55)
  colormap('default');  xlabel x, ylabel t, zlabel u
