grid_points = 5;
x = linspace(0, 1, 20 * grid_points);
xbar = linspace(0, 1, grid_points);

f = @(z) z;
f_true = @(z) (1/6) .* (z.^3) - (1/6) .* z;

% Plot Green's functions
%
clf;
y = zeros(grid_points, length(x));

for i = 1:1:grid_points
    y(i,:) = green_function(xbar(i), x, f(xbar(i)));
    plt = plot(x, y(i,:), '--', 'DisplayName', sprintf('G(x, %0.5f)', xbar(i)));
    hold on;
end

plt = plot(x, sum(y), 'DisplayName', 'F(x)');
plt = plot(x, f_true(x), 'DisplayName', 'u(x)');


title("Green's Functions");
xlabel('x');
ylabel('u(x)');
legend;

hold off;

waitfor(plt);
