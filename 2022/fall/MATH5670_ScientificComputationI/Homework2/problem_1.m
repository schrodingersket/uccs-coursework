#!/usr/bin/octave
a = 0;
b = 1;
grid_points = 5;
h = (b - a) / (grid_points - 1);
x = linspace(0, 1, 20 * grid_points);
xbar = linspace(0, 1, grid_points);

f = @(z) z;
f_true = @(z) ((1/6) .* (z.^3) - (1/6) .* z);

% Plot Green's functions
%
clf;
y = zeros(grid_points, length(x));

for i = 1:1:grid_points
    y(i,:) = h * f(xbar(i)) .* green_function(xbar(i), x);
    plt = plot(x, y(i,:), '--', 'DisplayName', sprintf('G(x, %0.2f)', xbar(i)));
    hold on;
end

plt = plot(x, sum(y), 'DisplayName', 'U(x)');
plt = plot(x, f_true(x), 'DisplayName', 'u(x)');


title("Green's Functions");
xlabel('x');
ylabel('u(x)');
legend;

hold off;

waitfor(plt);
