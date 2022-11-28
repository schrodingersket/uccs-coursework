grid_points = [50, 100, 200];
final_time = 0.1;
alpha = 0.4;

h = zeros(size(grid_points));
k = zeros(size(grid_points));
err = zeros(size(grid_points));

function eta = ic(x)
  % initial data

  beta = 600;
  eta = exp(-beta*(x - 0.5).^2);
  return
end

for n=1:length(grid_points)
    [h(n), k(n), err(n), _] = advection_lf_pbc(grid_points(n) - 1, alpha, final_time, @ic);
end

print('-dpng', 'problem_4a')

error_table(h, err);
error_loglog(h, err, 'problem_4a_error');