[X,Y] = meshgrid(linspace(0, 4, 50), linspace(0, 2, 50));

approximation_terms = 100;
u = @(x, y, n) (4/(n * pi * sinh(-2 * n * pi))) .* sinh(n * pi * (x - 4) / 2) .* sin(n * pi * y / 2);

figure;

Z = u(X, Y, 1);

if approximation_terms > 1
    for n=2:approximation_terms
        if mod(n, 2) ~= 0
            Z = Z + u(X, Y, n);
        end
    end
end

surf(X,Y,Z);
title('$\Delta u = 0$', 'interpreter', 'latex', 'VerticalAlignment', 'bottom');
colorbar();
view(2);
daspect([1 1])

print('-dpng', 'problem_1');

