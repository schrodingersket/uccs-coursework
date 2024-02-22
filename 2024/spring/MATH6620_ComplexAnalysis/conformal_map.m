% conformal_map.m
clear all;
clc;
clf;

hold off;
grid on;

t = linspace(0, 2*pi)
t1 = linspace(-3/2, 6);
t2 = linspace(-6, 6);
x = 2*t1+3;
y = -0.5*t2+3;
z1 = complex(t1, x);
z2 = complex(t2, y);
T_z1 = (z1 - i)./(z1 + i);
T_z2 = (z2 - i)./(z2 + i);
h_title = title('Map of orthogonal lines to unit disc with $T=\frac{z+i}{z-i}$', 'interpreter', 'latex');
set (h_title, 'fontsize', 16);
hold on;

plot(z1, 'DisplayName', '$z_1 = t + i(2t+3), \, t \in \left[-\frac{3}{2}, 6 \right]$')
plot(z2, 'DisplayName', '$z_2 = t + i\left(\frac{1}{2} t + 3\right), \, t \in [-6, 6]$')


xlim([-2*pi 2*pi]);
ylim([-1 2*pi]);
axis equal;


plot(T_z1, 'o', 'DisplayName', '$T(z_1)$')
plot(T_z2, 'o', 'DisplayName', '$T(z_2)$')
plot(exp(i*t), '--', 'DisplayName', 'Unit disc')

h_legend = legend('interpreter', 'latex');
set (h_legend, 'fontsize', 16);

print('MATH6620_ThompsonJ_ConformalMap', '-dpng')
