clc
clf

grid_points = [49];

x0 = 0;
x1 = 1;


% Boundary conditions (post-transform)
%
% beta = 1;
beta = airy(0, 0);

for i=1:length(grid_points)
    m = grid_points(i);
    fprintf('Computing solution for %d interior grid points...\n', m)
    
    h = (x1 - x0) / (m+1);

    % Create solution matrix A 
    %
    y_k = @(k) k*h;
    A = zeros(m+1, m+1);
    for n=1:m
        y_n = y_k(n);
        if n > 1
            A(n, n-1) = 1/h^2 - 1/(2*y_n*h);
        end
            
        A(n, n) = -2/h^2 + log(y_n) / y_n^2;

        if n < m + 1
            A(n, n+1) = 1/h^2 + 1/(2*y_n*h);
        end
    end
    A(m+1, m+1) = 1;

    % RHS
    %
    b = zeros(m, 1);
    b(m+1) = beta;

    % Solve system
    %
    u = A\b;

    % Plot (after denormalizing variables)
    %
    y = [x0:h:x1];
    y_int = y(2:length(y));
    x = -log(y_int);
    x_ref = linspace(0, max(x));
    u_true = airy(0, x_ref);

    hold on;
    plot(x, u, 'o', 'DisplayName', '$U(x)$');
    plot(x_ref, u_true, '--', 'DisplayName', '$Ai(x)$');
    title(sprintf('Airy function on $\\tilde{x} = [0, 1]$ with $h=%3.2f$', h), 'interpreter', 'latex')
    legend('interpreter', 'latex');

    print('-dpng', 'problem_4a_airy_function')
end
