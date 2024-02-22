step_sizes = [19, 29, 39, 49, 99];
h_vec = zeros(length(step_sizes), 1);
k_vec = zeros(length(step_sizes), 1);
err_vec = zeros(length(step_sizes), 1);


% true solution for comparison:
% For Gaussian initial conditions u(x,0) = exp(-beta * (x-0.4)^2)
beta = 150;
kappa = .02;               % heat conduction coefficient:
alpha = 4;                 % k = alpha * h
utrue = @(x,t) exp(-(x-0.4).^2 / (4*kappa*t + 1/beta)) / sqrt(4*beta*kappa*t+1);

for i=1:length(step_sizes)
    t_capture = NaN;
    if i == length(step_sizes)
        t_capture = -1;
    end
    [h_vec(i), k_vec(i), err_vec(i)] = heat_trbdf2(step_sizes(i), 0, 1, kappa, alpha, utrue, t_capture, 'problem_2b_');
end

disp('        h            k          error')
disp([h_vec, k_vec, err_vec]);

error_loglog(h_vec, err_vec, 'problem_2b_heatTRBDF2_error');