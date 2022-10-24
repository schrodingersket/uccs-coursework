clc
clf

steps = [5, 20, 100, 200, 300];

x0 = 0;
x1 = 1;

n_lambdas = 5;
lambdas = zeros(length(steps), n_lambdas);
step_sizes = (x1 - x0) ./ steps;

for i=1:length(steps)
    m = steps(i);
    fprintf('Calculating eigenvalues for %d interior points...\n', m)
    
    h = (x1 - x0) / m;

    % Create internal matrix A 
    %
    e = ones(m, 1);
    A = spdiags([e -2*e e], -1:1, m, m);

    alpha = 0;
    beta = 0;

    % Set boundary conditions
    %
    
    % Left Dirichlet condition
    %
    A(1, 1) = A(1, 1) + alpha*h^2;
    
    % Right Dirichlet condition
    %
    A(m, m) = A(m, m) + beta*h^2;

    % RHS
    B = spdiags(e, 0, m, m);
    %B(m, m) = 0;


    [eigvec, eigval] = eig(full(A), full(B));
    alpha = flip(diag(eigval));
    
    lambda = alpha ./ (h^2);
    lambdas(i, :) = lambda(1:n_lambdas);


    if i == length(steps)
        figure;
        title('Eigenvectors')
        
        for j = 1:n_lambdas
            plot(eigvec(:, j), 'DisplayName', sprintf('Eigenvector for j=%d', j));
            hold on;
        end
        
        hold off;
        legend;
        shg;
    end
end



disp(' ')
disp('       n     lambda_1 lambda_2  lambda_3  lambda_4  lambda_5')
disp([steps' lambdas])
fprintf('      %d      %3.4f   %3.4f   %3.4f  %3.4f  %3.4f\n', Inf, -(1*pi)^2, -(2*pi)^2, -(3*pi)^2, -(4*pi)^2, -(5*pi)^2);

eig_err = abs(lambdas(:, 1) - (-pi^2));
error_loglog(step_sizes, eig_err)

input('Press [Enter] to continue...')