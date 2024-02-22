grid_points = [50, 100, 100];
final_time = 0.1;
alpha = 1.0;

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
    [h(n), k(n), err(n)] = advection_sklf_pbc(grid_points(n) - 1, alpha, final_time, @ic);
end

print('-dpng', 'problem_1a')

error_table(h, err);
error_loglog(h, err, 'problem_1a_error');