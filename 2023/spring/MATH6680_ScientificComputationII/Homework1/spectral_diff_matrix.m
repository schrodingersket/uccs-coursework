function spectral_diff_matrix(ufn, uprimefn, filename)
    % p1.m - convergence of fourth-order finite differences
    % For various N, set up grid in [-pi,pi] and function u(x):
    Nvec = 2:2:100;

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
    
        % Construct spectral differentiation matrix:
        %
        tStart = tic;
        column = [0 .5*(-1).^(1:N-1) .* cot((1:N-1) * h/2)];
        D = toeplitz(column, column([1 N:-1:2]));

        % Plot max(abs(D*u-uprime)):
        %
        Nerror(i) = norm(D * u - uprime, inf);
        T(i) = toc(tStart);
    end
     
    grid on
    
    loglog(Nvec, Nerror, '.', 'markersize', 15, 'DisplayName', 'Spectral Differentiation')
    
    xlabel('N')
    ylabel('Error')
    title('Convergence of spectral differentiation')
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
