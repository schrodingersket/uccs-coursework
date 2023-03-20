#!/usr/bin/env octave

x_0 = 0;
x_1 = 1;
x_bar = linspace(x_0, x_1);
max_computation_order = 1;
u_true = @(x) x; % true initial condition
u_steady = (1/2) * ones(size(x_bar)); % true steady state

% Computes a particular Fourier series term
f_term = @(x, n, a_n, t) a_n * cos(n*pi*x) * exp(-t * (n*pi)^2);

% Computes a particular Fourier cosine coefficient
fs_coefficient = @(n) (1 - (-1)^(n)) * (-2) / (n*pi)^2;

for i=0:max_computation_order
    U = (1/2) * ones(size(x_bar));
    terms = 5 * 10^i;

    for j=1:terms
        U = U + f_term(x_bar, j, fs_coefficient(j), 0);
    end

    figure; 
    hold on;
    plot(x_bar, U, 'DisplayName', 'Fourier Series Approximation')
    plot(x_bar, u_true(x_bar), 'DisplayName', 'u(x, 0)')
    title(sprintf(
        'Fourier series approximation to u(x, 0) = x with %d terms',
         terms
    ))
    xlabel('x')
    ylabel('u(x, 0)')

    legend('Location', 'northwest') 
    print('-dpng', sprintf(
        'problem1_initial_fourier_series_solution_%d_terms',
         terms
    ))
    hold off;
    
    input('Press [Enter] to continue...');
end

for i=1:2
    U = (1/2) * ones(size(x_bar));
    terms = 2^(i^3) - 1;
    t = 0.2;

    for j=1:terms
        U = U + f_term(x_bar, j, fs_coefficient(j), t);
    end

    figure; 
    hold on;
    plot(x_bar, U, 'DisplayName', sprintf('u(x, %0.2f)', t))
    plot(x_bar, u_steady, 'DisplayName', 'u_s(x)')
    title(sprintf(
        'Fourier series approximation to u(x, %0.2f) with %d terms',
        t,
        terms
    ))
    xlabel('x')
    ylabel(sprintf('u(x, %0.2f)', t))
    axis([0 1 0 1])
    daspect([1 1 1])

    legend('Location', 'northwest') 
    print('-dpng', sprintf(
        'problem1_fourier_series_solution_%d_terms_t_%0.2f.png',
         terms,
         t
    ))
    hold off;
    
    input('Press [Enter] to continue...');
end