% clear all;

global epsilon k;

epsilon = 0.03;

N = 64;
dt = .4/N^2;
x = (2*pi/N)*(-N/2:N/2-1)';
k = [0:N/2-1 0 -N/2+1:-1]';

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


% RK4 integration
%
yn = V0;
tn = 0;

% Set to 1 to use integrating factor
%
integrating_factor = 1;

for n = 1:nmax
    if integrating_factor == 1
        k1 = fourier_burgers_if(tn, yn);
        k2 = fourier_burgers_if(tn + dt/2, yn + dt * k1/2);
        k3 = fourier_burgers_if(tn + dt/2, yn + dt * k2/2);
        k4 = fourier_burgers_if(tn + dt, yn + dt * k3);
    else 
        k1 = fourier_burgers(tn, yn);
        k2 = fourier_burgers(tn + dt/2, yn + dt * k1/2);
        k3 = fourier_burgers(tn + dt/2, yn + dt * k2/2);
        k4 = fourier_burgers(tn + dt, yn + dt * k3);
    end

    tn = n*dt;
    yn = yn + 1/6 * (k1 + 2*k2 + 2*k3 + k4) * dt;

    % Keep track of various time steps for plotting
    %
    if mod(n, nplt) == 0
        if integrating_factor == 1
            u_plot = [u_plot real(ifft(exp(-epsilon * k.^2 * (n-1)*dt) .* yn))];
        else
            u_plot = [u_plot real(ifft(yn))];
        end
        t_plot = [t_plot tn];
    end
end

% Plot results
%
figure;
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


%% RK4 RHS functions

function y = fourier_burgers(t, V)
    global epsilon k;

    % Burgers equation in Fourier space without integrating factor
    %
    y = (1i * k / 2) .* fft(real(ifft(V)).^2) - epsilon * k.^2 .* V;
end

function y = fourier_burgers_if(t, V)
    global epsilon k;

    % Burgers equation in Fourier space with integrating factor
    %
    y = (1i * k / 2) .* exp(epsilon .* k.^2 .* t) .* fft(real(ifft( exp(-epsilon * k.^2 * t) .* V)).^2);
end