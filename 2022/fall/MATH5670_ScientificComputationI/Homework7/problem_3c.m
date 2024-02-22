grid_points = 200;
final_time = 0.1;
alphas = [0.4, 0.2, 0.1];

h = zeros(size(grid_points));
k = zeros(size(grid_points));
err = zeros(size(grid_points));

for n=1:length(alphas)
    [h(n), k(n), err(n)] = advection_upwind_pbc(grid_points - 1, alphas(n), final_time);
    print('-dpng', sprintf('problem_3c_alpha-%.2f.png', alphas(n)))
end

