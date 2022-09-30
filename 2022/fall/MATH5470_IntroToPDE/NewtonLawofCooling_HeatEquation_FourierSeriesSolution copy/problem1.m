clf;
clc;

t = linspace(.01, 4*pi);

alpha = 1/4;
f = @(x) alpha ./ x;
g = @(x) tan(x);
y = tan(t);

figure;
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
