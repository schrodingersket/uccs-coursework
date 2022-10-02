clf;
clc;

t = linspace(.01, 4*pi);

alpha = 1/4;
f = @(x) alpha ./ x;
g = @(x) tan(x);
y = tan(t);

hold on;
grid on;

title('Transcendental equation solution for \alpha = 0.25');

% Plot alpha / lambda
plot(linspace(.01, 4*pi), f(linspace(.01, 4*pi)));

% Plot tan(lambda)
t1 = linspace(.01, pi/2 - .01);
plot(t1, g(t1), 'Color', 'red');
t2 = linspace(pi/2 + .01, 3*pi/2 - .01);
plot(t2, g(t2), 'Color', 'red');
t3 = linspace(3*pi/2 + .01, 5*pi/2 - .01);
plot(t3, g(t3), 'Color', 'red');
t4 = linspace(5*pi/2 + .01, 7*pi/2 - .01);
plot(t4, g(t4), 'Color', 'red');
t5 = linspace(7*pi/2 + .01, 9*pi/2 - .01);
plot(t5, g(t5), 'Color', 'red');

axis([0 4*pi -2 2]);
legend('\alpha / \lambda', 'tan(\lambda)');
xlabel('\lambda');

print('-dpng', 'problem1_transcendental_eigenvalues.png')


fun = @(lambda, x) (cos(lambda * x)).^2;
c_m = @(lambda) sin(lambda) / (lambda * integral(@(x) fun(lambda, x), 0, 1, 'RelTol',0,'AbsTol',1e-12));

% Found visually via Desmos: https://www.desmos.com/calculator/6drk08awls
%
lambda1 = .48;
lambda2 = 3.219;
lambda3 = 6.323;
lambda4 = 9.451;

fprintf('b_1 is (20 - T_0)(%0.5f) \n', c_m(lambda1));
fprintf('b_2 is (20 - T_0)(%0.5f) \n', c_m(lambda2));
fprintf('b_3 is (20 - T_0)(%0.5f) \n', c_m(lambda3));
fprintf('b_4 is (20 - T_0)(%0.5f) \n', c_m(lambda4));

T0 = 30;
u0 = 20;
x_0 = 0;
x_1 = 1;
x_bar = linspace(x_0, x_1);
max_computation_order = 1;
u_true = @(x) x; % true initial condition
u_steady = T0 * ones(size(x_bar)); % true steady state

% Computes a particular Fourier series term
f_term = @(x, lambda, a_n, t) a_n * cos(lambda*x) * exp(-t * (lambda)^2);

% Computes a particular Fourier cosine coefficient
lambdas = [lambda1, lambda2, lambda3, lambda4];
fs_coefficients = [(20 - T0)*c_m(lambda1), (20 - T0)*c_m(lambda2), (20 - T0)*c_m(lambda3), (20 - T0)*c_m(lambda4)];

terms = length(fs_coefficients);
    
figure; 
hold on;

title(sprintf('Fourier series approximation to u(x, t_0) with %d terms', terms));
xlabel('x');
ylabel('Temperature u(x, t_0)');

if u0 > T0
    axis([x_0 x_1 T0-0.5 u0+0.5]);
else
    axis([x_0 x_1 u0-0.5 T0+0.5]);
end

for t=0:2:10
    U = (T0) * ones(size(x_bar));
    for j=1:terms
        U = U + f_term(x_bar, lambdas(j), fs_coefficients(j), t);
    end

    plot(x_bar, U, 'DisplayName', sprintf('u(x, %0.2f)', t));
end

plot(x_bar, u_steady, 'DisplayName', 'u_s(x)');

legend('Location', 'southeast');
hold off;
print('-dpng', sprintf('problem1_fourier_series_solution_%d_terms_t0_%d.png', terms, T0));
