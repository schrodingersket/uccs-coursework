T = 1;

step_sizes = [19 29 39 49, 99];
h_vec = zeros(length(step_sizes), 1);
k_vec = zeros(length(step_sizes), 1);
err_vec = zeros(length(step_sizes), 1);

for i=1:length(step_sizes)
    [h_vec(i), k_vec(i), err_vec(i)] = heat_trbdf2(step_sizes(i));
end

disp('        h            k          error')
disp([h_vec, k_vec, err_vec]);

error_loglog(h_vec, err_vec, 'problem_2b_heatTRBDF2_error');