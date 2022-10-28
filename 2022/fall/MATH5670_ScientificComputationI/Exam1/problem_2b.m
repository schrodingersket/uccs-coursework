clc
clf

steps = [4, 10, 50, 100];
% steps = [4];

t0 = 1;
tf = 2;

u0 = 3;

u_true = @(t) t.^2 + 2.*t;
tt = linspace(t0, tf);

    
figure(1);
title("Numerical solution to $u'(t) = 2(t+1)$", 'interpreter', 'latex');

for i=1:length(steps)
    m = steps(i);
    k = (tf - t0) / (m+1);

    % Create internal matrix A 
    %
    e = ones(m, 1);
    A = 1/(4*k) * spdiags([-e 0*e e], -1:1, m, m);

    % RHS
    % 
    b = zeros(m, 1);
    b(1) = 2 + 1*k + u0/(4*k);
    for j=2:m-1
        b(j) = 2 + j*k;
    end
    b(m) = 2 + m*k + u0/(4*k);

    % disp(full(A));
    % disp(b);
    u = A\b;    

    subplot(2, 2, i);
    plot(linspace(t0 + k, tf-k, m), u, 'DisplayName', '$U(t)$');
    hold on;
    plot(tt, u_true(tt), 'DisplayName', '$u_{true}(t)$')
    title(sprintf('$\\Delta t = %3.2f$', k), 'interpreter', 'latex');
    legend('interpreter', 'latex', 'Location', 'southeast');
end

print('-dpng', 'problem_2b_solution.png')