clf;
clear all;
rng(0);  % fix random number generator to zero

t = (0:0.2:2*pi)';
p = (0:0.2:2*pi)';
[T, P] = meshgrid(t, p);
noiselevel = 0.1;

% Parametric Fernandez-Guasti squircle
%
r = 1;
s = 0.75;
X = r * cos(T) .* cos(P) ./ sqrt(1 - s * cos(T).^2 .* sin(P).^2 - s * sin(T).^2) + noiselevel * rand(size(t));
Y = r * cos(T) .* sin(P) ./ sqrt(1 - s * cos(T).^2 .* cos(P).^2 - s * sin(T).^2);
Z = r * sin(T) ./ sqrt(1 - s * cos(T).^2);

x = X(:);
y = Y(:);
z = Z(:);


N = length(T)*length(P); % total number of points in the point cloud

% Plot point cloud
%
figure(1), clf

plot3(x, y, z, 'k.')
title('3-D Point Cloud Source')

axis equal
grid on
hold on

print('problem_2ii_source_cloud', '-dpng')

xzero = [x y z];
xzero = [xzero(end, :); xzero; xzero(1, :)];

normal = zeros(N, 3);

% radius used to collect points for normal calculation
%
normal_radius = 0.30;
max_neighbor_radius = 0;
for j=1:size(normal, 1)
    data = xzero(j, :);
    min_radius = 3;
    for k = 1:size(xzero, 1)
        d = norm(xzero(j, :) - xzero(k, :));

        if d ~= 0 && d < min_radius
            min_radius = d;
        end

        if j ~= k && (0 < d) && (d < normal_radius)
            data = [data; xzero(k, :)];
        end
    end

    if min_radius > max_neighbor_radius
        max_neighbor_radius = min_radius;
    end

    normal(j,:) = bestfitnormal(data);
end

max_neighbor_radius

xzero([1, N + 2], :)=[];
alpha=4e-1;

% Exterior and interior points
%
xplus = xzero + alpha*normal;
xminus = xzero - alpha*normal;

% Plot interior/exterior point cloud
%
figure(2), clf
plot3(xzero(:, 1), xzero(:, 2), xzero(:, 3), 'ro', 'DisplayName', 'x_j')
hold on
plot3(xplus(:, 1), xplus(:, 2), xplus(:, 3), 'b*', 'DisplayName', 'x_j^+')
plot3(xminus(:, 1), xminus(:, 2), xminus(:, 3), 'gd', 'DisplayName', 'x_j^-')

axis equal
legend()
grid on


view([-90 0])
print('problem_2ii_normal_cloud_2d', '-dpng')
view([-25, 55])
print('problem_2ii_normal_cloud_3d', '-dpng')


eps = 10;
% phi = @(r) exp(-(eps*r).^2); %GA
% phi = @(r) 1./sqrt(1 + (eps*r).^2); %IMQ
phi = @(r) sqrt(1 + (eps*r).^2); %MQ

M = 3*N;
RHS = [zeros(N,1); alpha*ones(N,1); -alpha*ones(N,1)];

A = zeros(M);

Xbig = [xzero; xplus; xminus];

for j = 1:M
    for k = 1:M
        A(j,k) = phi(norm(Xbig(j,:)-Xbig(k,:)));
    end
end

C = A\RHS;

rbf = @(Xin) C' * phi(sqrt(sum((Xbig - repmat(Xin, size(Xbig, 1), 1)).^2, 2)));

[Xplot, Yplot, Zplot] = meshgrid(-1:0.2:1, -1:0.2:1, -1:0.2:1);

Uplot = zeros(size(Xplot));

for j = 1:size(Xplot, 1)
    for k = 1:size(Yplot, 1)
        for l = 1:size(Zplot, 1)
            Uplot(j, k, l)=rbf([Xplot(j, k, l) Yplot(j, k, l) Zplot(j, k, l)]);
        end
    end
end

figure(3), clf

isosurface(Xplot, Yplot, Zplot, Uplot, 0);
title('Level Surface F(x) = 0')
axis equal;
hold on;
grid on;
plot3(xzero(:, 1), xzero(:, 2), xzero(:, 3), 'k.')
view([-55, 25])
print('problem_2ii_rbf_surface', '-dpng')
