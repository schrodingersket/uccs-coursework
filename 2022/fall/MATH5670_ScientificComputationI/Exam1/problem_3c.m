
clc
clf

steps = [4, 10, 50, 100];

t0 = 1;
tf = 5;

u0 = 3;

% True solution, which in this case is known
%
u_true = @(t) t.^2 + 2.*t;
tt = linspace(t0, tf);

f = @(t) 2*(t+1);
    
figure(1);
title("Numerical solution to $u'(t) = 2(t+1)$", 'interpreter', 'latex');

for i=1:length(steps)
    m = steps(i);
    k = (tf - t0) / (m+1);
    t = linspace(t0, tf, m+1);
    tn = @(n) t0 + n*k;

    % Create internal matrix A 
    %
    u = u0 * ones(m+1, 1);

    for j=2:m+1
        Y1 = u(j-1) + (k/2) * f(tn(j-1) + k/2);
        u(j) = u(j-1) + k * f(tn(j-1) + k/2);
    end


    subplot(2, 2, i);
    plot(t, u, 'DisplayName', '$U(t)$');
    hold on;

    % Alternatively, we can plot against the true solution if known:
    %
    plot(tt, u_true(tt), '--', 'DisplayName', '$u_{true}(t)$')

    title(sprintf('$\\Delta t = %3.2f$', k), 'interpreter', 'latex');
    legend('interpreter', 'latex', 'Location', 'northwest');
end

print('-dpng', 'problem_3c_solution.png')