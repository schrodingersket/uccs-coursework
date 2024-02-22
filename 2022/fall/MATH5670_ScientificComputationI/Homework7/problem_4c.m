alpha = 0.4;
final_time = 0.1;

a = 2;           % advection velocity

ax = 0;
bx = 1;
m = 199;

h = (bx-ax)/(m+1);         % h = delta x
k = alpha*h;               % time step
nu = a*k/h;                % Courant number

global beta;
beta = 100;
global xi;
xi = 150;

function eta = ic(x)
  global beta;
  global xi;

  eta = exp(-beta*(x - 0.5).^2) .* sin(xi * x);
  return
end

group_velocity = a * cos(xi * h) / sqrt(1 - nu^2 * sin(xi * h)^2);

[_, _, _, w_avg] = advection_lf_pbc(m, alpha, final_time, @ic);

disp(sprintf('Predicted group velocity is %0.2f', group_velocity));
print('-dpng', 'problem_4c')