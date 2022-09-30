% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [a,b] x [a,b].  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)

function [hx, hy, err] = poisson_5pt(bx, ax, by, ay, m, n, f, u, plot_soln)
    arguments
        bx (1,1) double
        ax (1,1) double {mustBeLessThan(ax, bx)}
        by (1,1) double
        ay (1,1) double {mustBeLessThan(ay, by)}

        m (1,1) double {mustBePositive}
        n (1,1) double {mustBePositive}

        f  {mustBeA(f, "function_handle")}
        u  {mustBeA(u, "function_handle")}
        
        plot_soln (1, 1) logical
    end

    if ~exist('plot_soln','var')
        plot_soln = 0;
    end

    hx = (bx-ax)/(m+1);
    hy = (by-ay)/(n+1);

    % grid points for x including boundaries
    %
    x = linspace(ax, bx, m + 2);

    % Grid points for y including boundaries
    %
    y = linspace(ay, by, n + 2);
    
    
    [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
    X = X';                     % transpose so that X(i,j),Y(i,j) are
    Y = Y';                     % coordinates of (i,j) point
    
    Iint = 2:m+1;              % indices of interior points in x
    Jint = 2:n+1;              % indices of interior points in y
    Xint = X(Iint,Jint);        % interior points
    Yint = Y(Iint,Jint); 
    
    
    rhs = f(Xint,Yint);        % evaluate f at interior points for right hand side
                               % rhs is modified below for boundary conditions.

    % Evaluate true solution at grid points
    % 
    utrue = u(X, Y);
    
    % set boundary conditions around edges of usoln array:
    
    usoln = utrue;              % use true solution for this test problem
                                % This sets full array, but only boundary values
                                % are used below.  For a problem where utrue
                                % is not known, would have to set each edge of
                                % usoln to the desired Dirichlet boundary values.
    
    
    % adjust the rhs to include boundary terms:
    rhs(:,1) = rhs(:,1) - usoln(Iint, 1)/(hy^2);
    rhs(:,n) = rhs(:,n) - usoln(Iint, n + 2)/(hy^2);
    rhs(1,:) = rhs(1,:) - usoln(1, Jint)/(hx^2);
    rhs(m,:) = rhs(m,:) - usoln(m + 2, Jint)/(hx^2);
    
    
    % convert the 2d grid function rhs into a column vector for rhs of system:
    F = reshape(rhs, m*n, 1);
    
    % form matrix A:
    e = ones(max(m, n),1);
    T = spdiags([(1/hx^2)*e (-2/hx^2 - 2/hy^2)*e (1/hx^2)*e], [-1 0 1], m, m);
    S = spdiags([(1/hy^2)*e (1/hy^2)*e], [-1 1], n, n);
    A = (kron(speye(n), T) + kron(S, speye(m)));
    
    % Solve the linear system:
    uvec = A\F;  
    
    % reshape vector solution uvec as a grid function and 
    % insert this interior solution into usoln for plotting purposes:
    % (recall boundary conditions in usoln are already set) 
    
    usoln(Iint, Jint) = reshape(uvec, m, n);
    
    % assuming true solution is known and stored in utrue:
    err = max(max(abs(usoln - utrue)));   
    
    % plot results if specified:
    if plot_soln 
        clf
        hold on
        
        % plot grid:
        % plot(X,Y,'g');  plot(X',Y','g')
        
        % plot solution:
        contour(X, Y, usoln, 30, 'k')
        
        axis([ax bx ay by])
        daspect([1 1 1])
        title(sprintf('Contour plot of computed solution on [%0.2f, %0.2f] x [%0.2f, %0.2f]', ax, ay, bx, by))
        print('-dpng', sprintf('poisson_5pt_stencil_%0.0f-%0.0f_%0.0f-%0.0f_hx-%0.3f_hy-%0.3f.png', ax, ay, bx, by, hx, hy));
        hold off

        input('Press [Enter] to continue...')
    end
end