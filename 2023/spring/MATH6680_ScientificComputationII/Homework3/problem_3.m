% p.19.m - 2nd-order wave eq. on Chebyshev grid (compare p6.m)

% Time-stepping leap by leap frog formula:
N = 80
x = cos(pi*(0:N)/N);
dt = 8/N^2;
v = exp(-200*x.^2);
vold = exp(-200*(x-dt).^2);
tmax = 4;
tplot = 0.075;
plotgap = round(tplot/dt);
dt = tplot/plotgap;
nplots = round(tmax/tplot);
plotdata = [v; zeros(nplots, N + 1)];
tdata = 0;
clf
drawnow
h = waitbar(0, 'please wait...');

for i = 1:nplots, waitbar(i/nplots)
    for n = 1:plotgap
        w = chebfft(chebfft(v))';
        w(1) = 0;
        w(N+1) = 0;
        vnew = 2*v - vold + dt^2*w;
        vold = v;
        v = vnew;
    end
    plotdata(i+1, :) = v;
    tdata = [tdata; dt*i*plotgap];
end

% Plot results:
clf, drawnow, waterfall(x, tdata, plotdata)
axis([-1 1 0 tmax -2 2])
view(10, 70)
grid off
colormap([0 0 0])
ylabel t
zlabel u
close(h)

print('problem_3', '-dpng')