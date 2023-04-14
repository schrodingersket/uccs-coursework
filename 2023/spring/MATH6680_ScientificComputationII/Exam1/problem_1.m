clear all;

N = 64;
dt = .1/N^2;
x = (2*pi/N)*(-N/2:N/2-1)';

% Initial condition
%
u = exp(-10 * (sin(x/2)).^2);

% Initial condition in Fourier space (where we do the RK4 boogie)
%
V0 = fft(u);

K = [0:N/2-1 0 -N/2+1:-1]';

tmax = 1;
max_plot_lines = 25;
nplt = floor((tmax/max_plot_lines)/dt);
nmax = round(tmax/dt);

% Initial data for plotting
%
u_plot = u;
t_plot = 0;

% RHS of the ODE
%
function y = fourier_burgers(t, V)
    epsilon = 0.03;

    n = size(V, 1);
    k = [0:n/2-1 0 -n/2+1:-1]';

    y = (1i * k / 2) .* exp(epsilon .* k.^2 .* t) .* fft(real(ifft( exp(-epsilon * k.^2 * t) .* V)).^2);
end

% RK4 integration
%
yn = V0;
tn = 0;

for n = 1:nmax
    k1 = fourier_burgers(tn, yn);
    k2 = fourier_burgers(tn + dt/2, yn + dt * k1/2);
    k3 = fourier_burgers(tn + dt/2, yn + dt * k2/2);
    k4 = fourier_burgers(tn + dt, yn + dt * k3);

    tn = n*dt;
    yn = yn + 1/6 * (k1 + 2*k2 + 2*k3 + k4) * dt;

    % Keep track of various time steps for plotting
    %
    if mod(n, nplt) == 0
        u_plot = [u_plot real(ifft(yn))];
        t_plot = [t_plot tn];
    end
end

% Plot results
%
waterfall(x, t_plot, u_plot')
colormap([0 0 0])
view(-20, 25)
xlabel('x')
ylabel('t')
axis([-pi pi 0 tmax 0 1])
grid off
set(gca, 'ztick', [0, 1])
pbaspect([1 1 .13])

print('-dpng', 'problem_1')