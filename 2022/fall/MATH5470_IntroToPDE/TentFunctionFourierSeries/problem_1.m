x_0 = -1;
x_1 = 1;
L_0 = (x_1 - x_0)/2;
a_0 = 1/2;
t = linspace(x_0, x_1);
max_computation_order = 1;

u_true = @(x) 1 - abs(x);

% Computes a particular Fourier series term
f_term = @(x, n, L, a_n) a_n * cos((2*n - 1)*pi*x / L);

% Computes a particular Fourier cosine coefficient
fs_coefficient = @(n, L)  4 / ((2*n - 1)*pi/L)^2;

for i=0:max_computation_order
    U = a_0 * ones(size(t));
    terms = 1 * 5^i;
    
    for j=1:terms
        U = U + f_term(t, j, L_0, fs_coefficient(j, L_0));
    end

    figure; 
    hold on;
    plot(t, U, 'DisplayName', 'Fourier Series Approximation')
    plot(t, u_true(t), 'DisplayName', 'f(x)')
    title(sprintf('Fourier series approximation to f(x) = 1 - |x| with %d term(s)', terms))
    xlabel('x')
    ylabel('u_0(x)')

    legend('Location', 'northwest') 
    print('-dpng', sprintf('problem1_fourier_series_solution_%d_term', terms))
    hold off;
    
    input('Press [Enter] to continue...');
end
