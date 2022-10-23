steps = [5, 20, 30, 50];

x0 = 0;
x1 = 1;

n_lambdas = min(steps);
lambdas = zeros(length(steps), n_lambdas);
step_sizes = (x1 - x0) ./ steps;

for i=1:length(steps)
    m = steps(i);
    disp(sprintf('Calculating eigenvalues for %d interior points...', m))
    
    h = (x1 - x0) / m;

    % Create internal matrix A 
    %
    e = ones(m, 1);
    A = spdiags([e -2*e e], -1:1, m, m);

    alpha = 0;

    % Set boundary conditions
    %
    
    % Left Dirichlet condition
    %
    A(1, 1) = A(1, 1) + alpha;
    
    % Right Robin condition
    %
    A(m, m) = 4/(3 + 2*h) - 2;
    A(m, m-1) = 1 - 1/(3 + 2*h);

    [eigvec eigval] = eig(A);
    alpha = flip(diag(eigval));
    
    lambda = alpha ./ (h^2);
    lambdas(i, :) = lambda(1:n_lambdas);
    
    if i == length(steps)
        figure;
        hold on;
        for j = 1:n_lambdas
            plot(eigvec(:, j), 'DisplayName', sprintf('Eigenvector for j=%d', j));
        end

        legend;

        input('Press [Enter] to continue...')
        hold off;
    end
end



disp(' ')
disp('       n       lambda_1  lambda_2   lambda_3   lambda_4   lambda_5')
disp([steps' lambdas])
input('Press [Enter] to continue...');