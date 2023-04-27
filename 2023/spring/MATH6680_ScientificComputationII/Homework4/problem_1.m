clf
clear all

rng(0)

t = 0:0.1:2*pi;

N = length(t); % Number of points in point cloud
M = 2; % Number of points sampled to compute normal plane

x = 2*cos(t) + 0.1*rand(size(t));
x = [x(1) x x(end)];
y = 3*sin(t);
y = [y(1) y y(end)];

xzero = [x(2:end-1)' y(2:end-1)'];

hold on;
% plot(x, y, 'k.')

normal = zeros(N, 2);

for j = 2:N-1
    data = [ x(j-1:j+1)' y(j-1:j+1)' ];
    normal(j, :) = bestfitnormal(data);
end

axis equal

alpha = 0.1;

% Exterior/interior points
%
xplus = xzero + alpha*normal;
xminus = xzero - alpha*normal;

plot(xplus(:, 1), xplus(:, 2), 'b*')
plot(xminus(:, 1), xminus(:, 2), 'gd')

print('problem_1_point_cloud', '-dpng')