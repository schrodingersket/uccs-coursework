grid_points = [50, 100, 200];
final_time = 0.1;
alpha = 0.4;

h = zeros(size(grid_points));
k = zeros(size(grid_points));
err = zeros(size(grid_points));

for n=1:length(grid_points)
    [h(n), k(n), err(n)] = advection_upwind_pbc(grid_points(n) - 1, alpha, final_time);
end

print('-dpng', 'problem_3b')

error_table(h, err);
error_loglog(h, err, 'problem_3b_error');