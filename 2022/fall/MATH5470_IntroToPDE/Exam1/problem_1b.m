x_0 = 0;
x_1 = pi;
x_bar = linspace(x_0, x_1);
max_computation_order = 1;
u_true = @(x) 10 + x; % true initial condition
u_steady = @(x) 10; % true steady state

figure; 
hold on;
plot(x_bar, u_true(x_bar), 'DisplayName', 'u(x, 0)')
title('Initial Temperature Profile u(x, 0)')
xlabel('x')
ylabel('u(x, 0)')

legend('Location', 'northwest') 
print('-dpng', 'problem1b_initial_temperature_profile')
hold off;

input('Press [Enter] to continue...');