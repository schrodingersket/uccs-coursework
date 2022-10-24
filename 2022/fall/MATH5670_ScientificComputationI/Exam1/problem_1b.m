clc
clf

steps = [5, 50, 100, 200, 300, 400];

x0 = 0;
x1 = 1;

n_lambdas = min(steps);
lambdas = zeros(length(steps), n_lambdas);
step_sizes = (x1 - x0) ./ steps;

for i=1:length(steps)
    m = steps(i);
    fprintf('Calculating eigenvalues for %d points...\n', m)
    
    h = (x1 - x0) / (m+1);

    % Create solution matrix A 
    %
    e = ones(m+1, 1);
    A = 1/h^2*spdiags([e -2*e e], -1:1, m+1, m+1);

    alpha = 0;

    % Set boundary conditions
    %
    
    % Left Dirichlet condition
    %
    A(1, 1) = A(1, 1) + alpha * h^2;
    
    % Second-order Right Robin condition
    %
    A(m+1, m+1) = (4/(3 + 2*h) - 2);
    A(m+1, m) = (1 - 1/(3 + 2*h));


    % First-order Right Robin condition
    %
%     A(m+1,m+1) = (1+h)/h^2;
%     A(m+1,m) = -1/h^2;

    % RHS
    B = spdiags(e, 0, m+1, m+1);
    B(m+1, m+1) = 0;

    [unsorted_eigvec, unsorted_eigval] = eig(full(A), full(B));
    [d,ind] = sort(diag(unsorted_eigval));
    eigval = unsorted_eigval(ind,ind);
    eigvec = unsorted_eigvec(:,ind);


    alpha = flip(diag(eigval));
    
    lambda = alpha;
    lambdas(i, :) = lambda(1:n_lambdas);
    
    if i == length(steps)
        for j = 1:n_lambdas
            plot(eigvec(:, j), 'DisplayName', sprintf('Eigenvector for j=%d', j));
            title('Eigenvectors')
            hold on;
        end

        hold off;
        legend;
        shg;
    end
end



disp(' ')
disp('       n    lambda_1  lambda_2  lambda_3  lambda_4  lambda_5')
disp([steps' lambdas])
