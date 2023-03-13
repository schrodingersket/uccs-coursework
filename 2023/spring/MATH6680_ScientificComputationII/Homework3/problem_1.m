for N = 1:4
    % Generate your favorite grid point distribution and store in the vector 'x'
    %
    [x, _] = gauss(N);
    % [_, x] = cheb(N);

    % Placeholder matrices
    %
    a = zeros(1, N);
    D = zeros(N, N);

    % Generate coefficients and store in the vector 'a'
    %
    for k = 1:N
        a(k) = prod(x(k) - x([1:(k-1), (k+1):end]));
    end

    for m = 1:N
        % Off-diagonal entries
        %
        for n = 1:N
            D(m, n) = a(m) / (a(n) * (x(m) - x(n)));
        end

        % Diagonal entries
        %
        D(m, m) = sum(1 ./ (x(m) - x([1:(m-1), (m+1):end])));
    end

    fprintf('\n\nDifferentiation matrix for N = %d Legendre points: \n', N)
    disp(D);
end