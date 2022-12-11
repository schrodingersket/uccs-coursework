tt = linspace(0, 2*pi);
thetas = [0, 1/4, 1/2, 3/4, 1];
x_min = -5;
x_max = 5;
y_min = -5;
y_max = 5;

x = [x_min x_min x_max x_max];
y = [y_min y_max y_max y_min];

S = @(t, theta) (-1 + exp(i*t))./(1 - theta  + theta*exp(i.*t));
stable_color = [.7 .7 .7];
for j=1:length(thetas)

    hold on;
    if thetas(j) > 1/2
        fill(x, y, stable_color, 'DisplayName', 'Stable Region');
        fill(real(S(tt, thetas(j))), imag(S(tt, thetas(j))), 'w', 'DisplayName', 'Unstable Region');
    elseif thetas(j) == 1/2
        fill(x, y, 'w', 'DisplayName', 'Unstable Region');
        fill([x_min x_min 0 0], [y_min y_max y_max y_min], stable_color, 'DisplayName', 'Stable Region');
    else
        fill(x, y, 'w', 'DisplayName', 'Unstable Region');
        fill(real(S(tt, thetas(j))), imag(S(tt, thetas(j))), stable_color, 'DisplayName', 'Stable Region');
    end
    hold off;

    xlim([x_min, x_max])
    ylim([y_min, y_max])
    title(sprintf('\\partialS  for \\theta=%0.2f', thetas(j)));
    xlabel('Real Axis')
    ylabel('Imaginary Axis')
    legend;
    print('-dpng', sprintf('problem_4b_t_%0.2f.png', thetas(j)));
    clf;
end
