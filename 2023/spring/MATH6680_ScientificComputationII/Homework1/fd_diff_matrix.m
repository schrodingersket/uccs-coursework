function fd_diff_matrix(ufn, uprimefn, filename)
    % p1.m - convergence of fourth-order finite differences
    % For various N, set up grid in [-pi,pi] and function u(x):
    Nvec = 2.^(3:16);

    Nerror = zeros(size(Nvec));
    T = zeros(size(Nvec));
    
    clf;
    subplot('position', [.1 .4 .8 .5]);
    
    subplot(2, 1, 1);
    hold on;

    for i = 1:length(Nvec)
        N = Nvec(i);
        h = 2*pi/N;
        x = -pi + (1:N)' * h;
        u = ufn(x);
        uprime = uprimefn(x);
    
        % Construct sparse 4th-order differentiation matrix:
        %
        tStart = tic;
        e = ones(N, 1);
        D = sparse(1:N, [2:N 1], 2*e/3, N, N) - sparse(1:N, [3:N 1 2], e/12, N, N);
        D = (D - D')/h;
    
        % Store finite difference computation time and error
        %
        Nerror(i) = norm(D * u - uprime, inf);
        T(i) = toc(tStart);
    end
     
    grid on
    
    loglog(Nvec, Nerror, '.', 'markersize', 15, 'DisplayName', 'Finite Differences')
    semilogy(Nvec, Nvec.^(-4), '--', 'DisplayName', 'N^{-4}')
    
    xlabel('N')
    ylabel('Error')
    title('Convergence of 4th-order finite differences')
    legend;
    
    % Plot partitions vs. computation time
    %
    subplot(2, 1, 2)
    
    plot(Nvec, T)
    
    xlabel('N')
    ylabel('Time')
    
    legend;
    print('-dpng', filename)
end
