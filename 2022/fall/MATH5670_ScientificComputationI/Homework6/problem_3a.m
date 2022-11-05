T = 1;

step_sizes = [18, 38, 48, 78, 98];
h_vec = zeros(length(step_sizes), 1);
k_vec = zeros(length(step_sizes), 1);
err_vec = zeros(length(step_sizes), 1);


% true solution for comparison:
% For u0 = -heaviside(0) + 1/2 initial condition 
global kappa = .02;        % heat conduction coefficient:
alpha = 1;                 % k = alpha * h

function u = utrue(x, t)
    global kappa;

    u = (1/2)*erfc(x./sqrt(4 * kappa * t));
    u(isnan(u)) = 0;
end

for i=1:length(step_sizes)
    [h_vec(i), k_vec(i), err_vec(i)] = heat_CN(step_sizes(i), -1, 1, kappa, alpha, @utrue);
end

disp('        h            k          error')
disp([h_vec, k_vec, err_vec]);

error_loglog(h_vec, err_vec, 'problem_3a_heatCN_error');