% p19.m - 2nd-order wave eq. on Chebyshev grid (compare p6.m)
clear all
close all

% Time-stepping by leap frog formula:
N = 80; 
[D, xp] = cheb(N);
D2 = D^2;

% Neumann condition u_x(-1, t) = 0
%
D2(N+1, :) = D(N+1, :);

x = xp';
dt = 8/N^2;
v = 0 * ones(size(x));
vold = 0 * ones(size(x - dt));
tmax = 5;
tplot = .025;
plotgap = round(tplot/dt); 
dt = tplot/plotgap;
nplots = round(tmax/tplot);
plotdata = [v; zeros(nplots,N+1)]; 
tdata = 0;
clf
drawnow
for i = 1:nplots
    for n = 1:plotgap
        w = (D2 * v')';
        w(1) = 0;
        w(N+1) = 0;
        vnew = 2*v - vold + dt^2*w;
        vold = v;
        v = vnew;

        % Dirichlet condition
        %
        t = dt*i*plotgap;
        v(1) = sin(10*t);

        % Neumann condition
        %
        v(N+1) = -1 / D(N+1, N+1) * (D(N+1, 1:N) * v(1:N)');
    end
    plotdata(i+1,:) = v;
    tdata = [tdata; dt*i*plotgap];
end

% Plot results:
clf
drawnow
waterfall(x,tdata,plotdata)
axis([-1 1 0 tmax -3 3])
view(10, 70)
grid off
colormap(1e-6*[1 1 1]); 
ylabel t
zlabel u

print('-dpng', 'problem_3.png')
pause
