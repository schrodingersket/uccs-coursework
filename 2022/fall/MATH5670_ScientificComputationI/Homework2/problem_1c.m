#!/usr/bin/octave
a = 0;
b = 1;
grid_points = 5;
h = (b - a) / (grid_points - 1);
x = linspace(a, b, 20 * grid_points);
xbar = linspace(a, b, grid_points);
alpha = 0;
beta = 0;

f = @(z) z;
f_true = @(z) ((1/6) .* (z.^3) - (1/6) .* z);

% Plot Green's functions
%
clf;
y = zeros(grid_points, length(x));

for i = 1:1:grid_points
    y(i,:) = h * f(xbar(i)) .* dirichlet_green_function(xbar(i), x);
    plt = plot(x, y(i,:), '--', 'DisplayName', sprintf('h * f(%0.2f) * G(x, %0.2f)', xbar(i), xbar(i)));
    hold on;
end

% Plot discrete approximation U
%
plt = plot(x, alpha*(1 - x) + beta*x + sum(y), 'DisplayName', 'U(x)');

% Plot exact solution for alpha = beta = 0
%
if alpha == 0 && beta == 0
    plt = plot(x, f_true(x), 'DisplayName', 'u(x)');
end


title("Green's Functions");
xlabel('x');
ylabel('u(x)');
legend;

hold off;

print('-dpng', 'problem1c_discrete_approx.png')
waitfor(plt);
