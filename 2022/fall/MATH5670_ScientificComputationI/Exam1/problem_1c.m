clf
clc

n_lambdas = 2;
lambdas = zeros(n_lambdas, n_lambdas*n_lambdas);
steps = [10, 15, 20, 25, 30];

a = 0;
b = 1;

step_sizes = (b - a) ./ steps;

function [lambda, eigvec] = poisson_5pt_eigval(bx, ax, by, ay, m, n)
    arguments
        bx (1,1) double
        ax (1,1) double {mustBeLessThan(ax, bx)}
        by (1,1) double
        ay (1,1) double {mustBeLessThan(ay, by)}

        m (1,1) double {mustBePositive}
        n (1,1) double {mustBePositive}
    end

    hx = (bx-ax)/(m+1);
    hy = (by-ay)/(n+1);
    n_lambdas = 4;
    
    % form matrix A:
    e = ones(max(m, n),1);
    T = spdiags([(1/hx^2)*e (-2/hx^2 - 2/hy^2)*e (1/hx^2)*e], [-1 0 1], m, m);
    S = spdiags([(1/hy^2)*e (1/hy^2)*e], [-1 1], n, n);
    A = (kron(speye(n), T) + kron(S, speye(m)));
    
    % Calculate eigenvalues:
    [unsorted_eigvec, unsorted_eigval] = eig(full(A));
    [d,ind] = sort(diag(unsorted_eigval));
    eigval = unsorted_eigval(ind,ind);
    eigvec = unsorted_eigvec(:,ind);

    lambda = flip(diag(eigval));
    
    % reshape eigenvector as a grid function for plotting purposes:
    % Iint = 1:m;
    % Jint = 1:n;
    % usoln(Iint, Jint) = reshape(eigvec, m, n);

    % assuming true solution is known and stored in utrue:
    % err = max(max(abs(usoln - utrue)));   
end

for k=1:length(steps)
    M = steps(k);
    N = steps(k);
    [eigenvalues eigenvectors] = poisson_5pt_eigval(b, a, b, a, M, N);
    lambdas(k, :) = eigenvalues(1:(n_lambdas*n_lambdas));
    
    % plot results if specified:
    if k == length(steps)
        hold on

        % grid points for x including boundaries
        %
        x = linspace(a, b, M);
    
        % Grid points for y including boundaries
        %
        y = linspace(a, b, N);
        
        [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
        X = X';                     % transpose so that X(i,j),Y(i,j) are
        Y = Y';                     % coordinates of (i,j) point
        
        % Plot eigenvectors
        %
        for l=1:n_lambdas^2
            subplot(2, 2, l)
            surf(X, Y, reshape(eigenvectors(:, l), M, N))
            title(sprintf('Eigenvector for n=%d', l))
        end

        print('-dpng', 'problem_1c_eigenvectors')
    end
end

disp(' ')
disp('      n    lambda_1 lambda_2 lambda_3 lambda_4')
disp([steps' lambdas])
fprintf('    %d    %3.3f  %3.3f  %3.3f  %3.3f\n', Inf, -(1^2 + 1^2)*pi^2, -(1^2 + 2^2)*pi^2, -(2^2 + 1^2)*pi^2, -(2^2 + 2^2)*pi^2);
disp(' ')

eig_err = abs(lambdas(:, 1) - (-(1 + 1)*pi^2));
error_loglog(step_sizes, eig_err)

print('-dpng', 'problem_1c_error.png')
