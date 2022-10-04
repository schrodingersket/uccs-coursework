#!/usr/bin/env octave

x_0 = 0;
x_1 = pi;
x_bar = linspace(x_0, x_1);
max_computation_order = 1;
u_true = @(x) 10 + x; % true initial condition
u_steady = 10 * ones(size(x_bar)); % true steady state

lambda = @(k, L) (2*k - 1) * pi / (2*L);

% Computes a particular Fourier series term
f_term = @(x, n, a_n, t) a_n * sin(lambda(n, x_1) * x) * exp(-t * (lambda(n, x_1))^2);

% Computes a particular Fourier cosine coefficient
fs_coefficient = @(n, L) (8 * L) / ((2*n - 1)^2 * pi^2) * (-1)^(n+1);


printf('First Fourier coefficient: %0.4f\n', fs_coefficient(1, x_1));
printf('Second Fourier coefficient: %0.4f\n', fs_coefficient(2, x_1));
printf('Third Fourier coefficient: %0.4f\n', fs_coefficient(3, x_1));
printf('Third Fourier coefficient: %0.4f\n', fs_coefficient(4, x_1));


% Plot Fourier series approximation to initial condition u(x, 0)
%
figure; 
hold on;
plot(x_bar, u_true(x_bar), 'DisplayName', 'u(x, 0)')
title('Fourier series approximation to u(x, 0) = 10 + x');
xlabel('x')
ylabel('u(x, 0)')

for i=0:max_computation_order
    U = 10 * ones(size(x_bar));
    terms = 1 * 10^i;

    for j=1:terms
        U = U + f_term(x_bar, j, fs_coefficient(j, x_1), 0);
    end

    plot(x_bar, U, 'DisplayName', sprintf('u_{%d}(x, 0)', terms));
end


legend('Location', 'northwest') 
print('-dpng', 'problem1e_initial_fourier_series')
hold off;

% Plot u(x, t)
%
for t=2:2:8
    figure; 
    hold on;
    title(sprintf('Fourier series approximation at t=%d', t))
    xlabel('x')
    ylabel(sprintf('u(x, %d)', t))
    axis([x_0 x_1 0 15])
    
    
    plot(x_bar, u_steady, 'DisplayName', 'u_s(x)')
    
    for i=0:1
        U = 10 * ones(size(x_bar));
        terms = 1 * 10^i;
    
        for j=1:terms
            U = U + f_term(x_bar, j, fs_coefficient(j, x_1), t);
        end
    
        plot(x_bar, U, 'DisplayName', sprintf('u_{%d}(x, %d)', terms, t))
    end
    
    legend('Location', 'northwest') 
    print('-dpng', sprintf( 'problem1e_fourier_series_solution_t_%d.png', t))
    hold off;
end

input('Press [Enter] to continue...');