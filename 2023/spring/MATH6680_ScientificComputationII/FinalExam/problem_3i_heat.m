clf;
clear all;

% Gaussian RBF
%
rbf = @(e, r) exp(-(e*r).^2);
d2rbf = @(e, r) 2 * e^2 * (2 * (e*r).^2 - 1) .* rbf(e, r);

% Thermal diffusivity
%
c = 0.5;

% Initial condition amplitude
%
T = 1;

% Max time
%
tmax = 1;

% Time plot interval
%
tplot = 0.1;

N = 20;

[D2, x] = D2RBF(N, rbf, d2rbf);

% Time step with explicit Euler formula
%
dt = min([10e-4, N^(-4)/c]); 
t = 0;
nplots = round(tmax / tplot);
plotgap = round(tplot / dt);
dt = tplot / plotgap;

% Initial condition (enforces hard initial BC)
%
v = T * cos(pi*x/2);

xx = -1:0.025:1;
vv = polyval(polyfit(x, v, N), xx);
plotdata = [vv; zeros(nplots, length(xx))];
tdata = t;


% Euler time step
for i = 1:nplots
    for n = 1:plotgap
        t = t + dt;
        v = v + dt * (c * (D2*v));  % Euler step for heat equation
        v(1) = 0;
        v(end) = 0;
    end

    vv = polyval(polyfit(x, v, N), xx);
    plotdata(i+1, :) = vv;
    tdata = [tdata; t];
end

% Generate plot
%
surf(xx, tdata, plotdata);
view(-35, 55);
colormap('default');
xlabel('x')
ylabel('t');
zlabel('u');

print('problem_3i_helmholtz_3d', '-dpng')

view(0, 90)

print('problem_3i_helmholtz_2d', '-dpng')
