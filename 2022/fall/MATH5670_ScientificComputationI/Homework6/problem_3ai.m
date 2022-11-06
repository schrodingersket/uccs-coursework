m = 39;

% true solution for comparison:
% For u0 = -heaviside(0) + 1/2 initial condition 
global kappa = .02;        % heat conduction coefficient:
alpha = 4;                 % k = alpha * h

function u = utrue(x, t)
    global kappa;

    u = (1/2)*erfc(x./sqrt(4 * kappa * t));
    u(isnan(u)) = 0;
end

heat_CN(m, -1, 1, kappa, alpha, @utrue, 1, 'problem_3ai_');
heat_CN(m, -1, 1, kappa, alpha, @utrue, -1, 'problem_3ai_');