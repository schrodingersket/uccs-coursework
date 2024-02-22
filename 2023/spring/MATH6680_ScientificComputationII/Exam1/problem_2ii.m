% Fisher equation u_t = u_xx + u(1 - u) on large domain [-L, L] 
% with front wave solution.
%
clear all;

function [xx, tdata, plotdata] = fisher_pde(L)
    clf;
    N = 32;
    [D, x] = cheb(N);
    D2 = D^2;
    D2([1 N+1], :) = zeros(2, N+1);
    
    dt = .01*N^(-2);
    dx = 1/N;
    
    t = 0;
    v = (1 + exp((L*x)/sqrt(6))).^(-2);
    
    tmax = 10;
    tplot = 0.5;
    nplots = round(tmax/tplot);
    plotgap = round(tplot/dt);
    dt = tplot/plotgap;
    
    xx = -1:dx:1;
    vv = polyval(polyfit(x, v, N), xx);
    
    plotdata = [vv; zeros(nplots, length(xx))];
    tdata = t;
    
    for i = 1:nplots
        for n = 1:plotgap
            t = t + dt;
            v = v + dt*((1/L)^2*D2*v + v - v.^2);
            
            % Enforce Dirichlet BCs explicitly
            %
            v(1) = 0;
            v(end) = 1;
        end
    
        plotdata(i + 1, :) = polyval(polyfit(L*x,v,N),L*xx);
        tdata = [tdata; t];
    end
    
    waterfall(L*xx, tdata, plotdata);
    grid on;
    % axis([-L L 0 tmax -0.5 1.5]);
    view(20, 25);
    colormap([0 0 0])
    xlabel('x')
    ylabel('t')
    zlabel('u')
    
    print('-dpng', sprintf('problem_2ii_L%d', L))
end

fisher_pde(20);
fisher_pde(40);
fisher_pde(60);