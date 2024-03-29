clc
clf

steps = [5, 50, 100, 200, 300, 400];

x0 = 0;
x1 = 1;

n_lambdas = 4;
lambdas = zeros(length(steps), n_lambdas);

for i=1:length(steps)
    m = steps(i);
    h = (x1 - x0) / (m+1);

    % Create internal matrix A 
    %
    e = ones(m, 1);
    A = 1/h^2*spdiags([e -2*e e], -1:1, m, m);

    % Set boundary conditions
    %
    alpha = 0;
    
    % Left Dirichlet condition
    %
    A(1, 1) = A(1, 1) + h^2*alpha;
    
    % Right Robin condition (second-order approx)
    %
    A(m, m) = (1/h^2) * (-2 + 4/(3 + 2*h));
    A(m, m-1) = (1/h^2) * (1 - 1/(3 + 2*h));

    [unsorted_eigvec, unsorted_eigval] = eig(full(A));
    [d,ind] = sort(diag(unsorted_eigval));
    eigval = unsorted_eigval(ind,ind);
    eigvec = unsorted_eigvec(:,ind);

    lambda = flip(diag(eigval));
    lambdas(i, :) = lambda(1:n_lambdas);

    if i == length(steps)
        figure(1);
        title('Eigenvectors');

        for j = 1:n_lambdas
            subplot(2, 2, j);
            plot(linspace(0, 1, length(eigvec(:, j))), eigvec(:, j));
            title(sprintf('Eigenvector for j=%d', j));
        end
        
        print('-dpng', 'problem_1b_eigenvectors.png')
    end
end


disp(' ')
disp('       n       lambda_1  lambda_2   lambda_3  lambda_4')
disp([steps' lambdas])