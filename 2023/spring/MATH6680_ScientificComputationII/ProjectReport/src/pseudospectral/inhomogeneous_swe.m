% Solve  the system [h ]   + [uh           ]   = [0     ]
%                   [uh]_t + [hu^2 + g'h^2/2]_x = [-g'hB_x]
% where g' = g cos(alpha)
%
clear all;
N = 2^8; % Number of spatial points
g = 9.81; % gravity
L = 10; % Interval length
b = L / (2*pi); % rescaling factor

t0 = 0; % initial time
tFinal = 4; % final time

dt = 1e-3; % time step length
nT = floor(tFinal/dt); % maximum number of time steps

x = linspace(0,L,N+1); % interval
x(end) = []; % Remove final point for periodic boundary conditions

h = 1/2 + 2/5 * sin(pi*x/L);

% Bathymetry function
%
B = @(x,t) 1/4*cos(2*pi*x/L) + 1/4;
Bx = @(x,t) -1/2 * pi/L * sin(2*pi*x/L);

alpha = @(x,t) atan(Bx(x, t));
gprime = @(x,t) g*cos(alpha(x,t));
% gprime = @(x,t) g;

u = zeros(size(h)); %initial velocity
uh = u.*h; %initial momentum/discharge (?)
h0 = h;
k = [0:N/2 - 1 0 -N/2+1:-1]; %fourier domain


RHS_h = @(h,uh,t) -real(ifft(1i*k.*fft(uh)));
RHS_uh = @(h, uh,t) -real(ifft(1i*k.*fft(uh.^2 ./ h + 1/2 * gprime(x, t) .* h.^2))) - gprime(x, t) .* h .* Bx(x,t);

tRan = linspace(0, tFinal, nT+1);

plotdatah = zeros(nT + 1, N);
plotdata_uh = zeros(nT + 1, N);
plotdata_u = zeros(nT + 1, N);

plotdata_uh(1,:) = uh;
plotdatah(1,:) = h;
u = uh./h;
plotdata_u(1,:) = u;

for jj = 1:nT
    t = jj*dt;
    k1 = (dt/b)*RHS_h(h, uh,t);
    k2 = (dt/b)*RHS_h(h + k1/2,uh,t);
    k3 = (dt/b)*RHS_h(h + k2/2,uh,t);
    k4 = (dt/b)*RHS_h(h + k3, uh,t);

    h = h + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

    k1 = (dt/b)*RHS_uh(h,uh,t);
    k2 = (dt/b)*RHS_uh(h,uh + k1/2,t);
    k3 = (dt/b)*RHS_uh(h, uh + k2/2,t);
    k4 = (dt/b)*RHS_uh(h, uh + k3,t);

    uh = uh + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    plotdatah(jj+1,:) = h;
    plotdata_uh(jj+1,:) = uh;
    u = uh./h;

    plotdata_u(jj + 1,:) = u;
end

% Plot results; since solutions are saved as CSV files, we leave these commented since we assume
% that the consuming script will likely plot these results differently regardless.
%
% subplot(1,2,1)
% waterfall(x, tRan, plotdatah)
% xlabel('$x$', 'fontsize', 25, 'interpreter', 'latex')
% ylabel('$t$', 'fontsize', 25, 'interpreter', 'latex')
% title('$h(t)$', 'fontsize', 20, 'interpreter', 'latex')
% hold on

% Uncomment these lines to also plot bathymetry function
%
% Bplot = B(x,t).'*ones(1,length(tRan));
% waterfall(x,tRan,Bplot.')
% colorbar

% subplot(1,2,2)
% waterfall(x, tRan, plotdata_u)
% xlabel('$x$', 'fontsize', 25, 'interpreter', 'latex')
% ylabel('$t$', 'fontsize', 25, 'interpreter', 'latex')
% title('$u(t)$', 'fontsize', 20, 'interpreter', 'latex')
% colorbar

[X,T] = meshgrid(x,tRan);
[AL, ~] = meshgrid(alpha(x,t), tRan);
ALcol = reshape(AL, length(x)*length(tRan),1);

Xcol = reshape(X, length(x)*length(tRan), 1);
Tcol = reshape(T, length(tRan)*length(x), 1);

uCol = reshape(plotdata_u, length(x)*length(tRan), 1);
hCol = reshape(plotdatah, length(x)*length(tRan), 1);

% Write CSV file
%
filename = '../data/inhomogeneous_swe_observation.csv';
fid = fopen(filename, 'wt');
fprintf(fid, '%s,%s,%s,%s,%s\n', 'x','t','v', 'h', 'alpha');  % header
fclose(fid);

dlmwrite('../data/inhomogeneous_swe_observation.csv', [ Xcol Tcol uCol hCol ALcol ], '-append');

