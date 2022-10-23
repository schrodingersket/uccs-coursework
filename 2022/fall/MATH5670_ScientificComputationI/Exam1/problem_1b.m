% steps = [5, 20, 30, 50];
steps = [2:20];

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
    A(m, m) = (-2 - 4*h)/(3 + 2*h);
    A(m, m-1) = (2 + 2*h)/(3 + 2*h);

    [eigval eigvec] = eig(A);
    alpha = flip(eigval);
    
    lambda = alpha ./ (h^2);
    lambdas(i, :) = lambda(1:n_lambdas);
    full(A)
end



disp(' ')
disp('       n       lambda_1  lambda_2   lambda_3   lambda_4   lambda_5')
disp([steps' lambdas])
input('Press [Enter] to continue...');