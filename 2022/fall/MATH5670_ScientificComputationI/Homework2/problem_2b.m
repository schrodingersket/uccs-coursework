#!/usr/bin/octave
a = 0;
b = 1;
grid_points = 5;
h = (b - a) / (grid_points - 1);
x = linspace(a + h, b - h, grid_points - 2);
xbar = linspace(a, b, grid_points);

% Create anonymous functions for Green's functions for boundary terms
%
G_0 = @(x) x - 1;
G_1 = @(x) ones(size(x));

B_interior = zeros(grid_points, length(x));

% Assign interior values
%
for i = 1:1:grid_points
    B_interior(i,:) = h * neumann_green_function(xbar(i), x);
end

% Print results
%
B = [ G_0(xbar)' B_interior G_1(xbar)' ];
display(B)

A = inv(B);
display(A)
