clc
clf

steps = [5, 20, 100, 200, 300];

x0 = 0;
x1 = 1;

n_lambdas = 4;
lambdas = zeros(length(steps), n_lambdas);
step_sizes = (x1 - x0) ./ steps;

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
    beta = 0;
    
    % Left Dirichlet condition
    %
    A(1, 1) = A(1, 1) + h^2*alpha;
    
    % Right Dirichlet condition
    %
    A(m, m) = A(m, m) + h^2*beta;

    [unsorted_eigvec, unsorted_eigval] = eig(full(A));
    [d,ind] = sort(diag(unsorted_eigval));
    eigval = unsorted_eigval(ind,ind);
    eigvec = unsorted_eigvec(:,ind);

    lambda = flip(diag(eigval));
    lambdas(i, :) = lambda(1:n_lambdas);

    if i == length(steps)
        figure(1);
        title('Eigenvectors')

        for j = 1:n_lambdas
            subplot(2, 2, j, 'align');
            plot(linspace(0, 1, length(eigvec(:, j))), eigvec(:, j));
            title(sprintf('Eigenvector for j=%d', j));
        end

        legend;
        shg;

        print('-dpng', 'problem_1a_eigenvectors.png')
    end
end



disp(' ')
disp('       n       lambda_1  lambda_2   lambda_3  lambda_4')
disp([steps' lambdas])
fprintf('      %d      %3.4f   %3.4f   %3.4f  %3.4f  %3.4f\n', Inf, -(1*pi)^2, -(2*pi)^2, -(3*pi)^2, -(4*pi)^2);

eig_err = abs(lambdas(:, 1) - (-pi^2));
error_loglog(step_sizes, eig_err)

print('-dpng', 'problem_1a_error.png')