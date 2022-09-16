x_0 = 0;
x_1 = pi;
t = linspace(x_0, x_1);
max_computation_order = 2;

u_true = @(x) x;

% Computes a particular Fourier series term
f_term = @(x, n, L, a_n) a_n * sin(n*pi*x / L);

% Computes a particular Fourier sine coefficient
fs_coefficient = @(n, L)  (-1)^(n+1) * (2*L) / (n*pi);

for i=1:max_computation_order
    U = zeros(size(t));
    terms = 5 * 10^i;
    
    for j=1:terms
        U = U + f_term(t, j, x_1, fs_coefficient(j, x_1));
    end

    figure; 
    hold on;
    plot(t, U, 'DisplayName', 'Fourier Series Approximation')
    plot(t, u_true(t), 'DisplayName', 'u_0(x)')
    title(sprintf(
        'Fourier series approximation to u_0(x) = x with %d terms',
         terms
    ))
    xlabel('x')
    ylabel('u_0(x)')

    legend('Location', 'northwest') 
    print('-dpng', sprintf(
        'problem1_fourier_series_solution_%d_terms',
         terms
    ))
    hold off;
    
    input('Press [Enter] to continue...');
end
