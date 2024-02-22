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

    
figure(1);
title("Numerical solution to $u'(t) = 2(t+1)$", 'interpreter', 'latex');

for i=1:length(steps)
    m = steps(i);
    k = (tf - t0) / (m+1);
    t = linspace(t0, tf, m+2);

    % Solve numerically to get required second point
    %
    [t_45, u_45] = ode45(@(t, u) 2*(t+1), t, u0);
    u1 = u_45(2);

    % Alternatively, we could find our needed second point via the analytical solution if known:
    %
    % u1 = u_true(t(2));

    % Create internal matrix A 
    %
    e = ones(m+2, 1);
    A = 1/(4*k) * spdiags([-e 0*e e], -2:0, m+2, m+2);

    % RHS
    % 
    b = zeros(m+2, 1);
    b(1) = u0/(4*k);
    b(2) = u1/(4*k);
    for j=3:m+2
        b(j) = 2 + j*k;
    end

    % disp([full(A) b]);

    u = A\b;    

    subplot(2, 2, i);
    plot(t, u, 'DisplayName', '$U(t)$');
    hold on;
    % Plot against ode45 solution
    %
    plot(t_45, u_45, 'o', 'DisplayName', '$u_{ode23}(t)$')

    % Alternatively, we can plot against the true solution if known:
    %
    plot(tt, u_true(tt), '--', 'DisplayName', '$u_{true}(t)$')

    title(sprintf('$\\Delta t = %3.2f$', k), 'interpreter', 'latex');
    legend('interpreter', 'latex', 'Location', 'northwest');
end

print('-dpng', 'problem_2b_solution.png')