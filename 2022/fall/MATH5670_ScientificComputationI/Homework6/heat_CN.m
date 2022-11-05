function [h,k,error] = heat_CN(m, ax, bx, kappa, alpha, utrue)
  %
  % heat_CN.m
  %
  % Solve u_t = kappa * u_{xx} on [ax,bx] with Dirichlet boundary conditions,
  % using the Crank-Nicolson method with m interior points.
  %
  % Returns k, h, and the max-norm of the error.
  % This routine can be embedded in a loop on m to test the accuracy,
  % perhaps with calls to error_table and/or error_loglog.
  %
  % From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

  clf              % clear graphics
  % hold on          % Put all plots on the same graph (comment out if desired)

  tfinal = 1;                % final time

  h = (bx-ax)/(m+1);         % h = delta x
  x = linspace(ax,bx,m+2)';  % note x(1)=0 and x(m+2)=1
  % u(1)=g0 and u(m+2)=g1 are known from BC's
  k = alpha*h;                  % time step

  nsteps = round(tfinal / k);    % number of time steps
  nplot = 1;      % plot solution every nplot time steps
  % (set nplot=2 to plot every 2 time steps, etc.)
  % nplot = nsteps;  % only plot at final time

  if abs(k*nsteps - tfinal) > 1e-5
    % The last step won't go exactly to tfinal.
    disp(' ')
    disp(sprintf('WARNING *** k does not divide tfinal, k = %9.5e',k))
    disp(' ')
  end


  % initial conditions:
  u0 = utrue(x,0);


  % Each time step we solve MOL system U' = AU + g using the Trapezoidal method

  % set up matrices:
  r = (1/2) * kappa* k/(h^2);
  e = ones(m,1);
  A = spdiags([e -2*e e], [-1 0 1], m, m);
  A1 = eye(m) - r * A;
  A2 = eye(m) + r * A;


  % initial data on fine grid for plotting:
  xfine = linspace(ax,bx,1001);
  ufine = utrue(xfine,0);

  % initialize u and plot:
  tn = 0;
  u = u0;

  plot(x,u,'b.-', xfine,ufine,'r')
  legend('computed','true')
  title('Initial data at time = 0')
  axis([ax bx 0 1])

  % input('Hit <return> to continue  ');


  % main time-stepping loop:

  for n = 1:nsteps
    t_next = tn + k;   % = t_{n+1}

    % boundary values u(0,t) and u(1,t) at times tn and t_next:

    g0n = u(1);
    g1n = u(m+2);
    g0n_next = utrue(ax,t_next);
    g1n_next = utrue(bx,t_next);

    % compute right hand side for linear system:
    uint = u(2:(m+1));   % interior points (unknowns)
    rhs = A2*uint;
    % fix-up right hand side using BC's (i.e. add vector g to A2*uint)
    rhs(1) = rhs(1) + r*(g0n + g0n_next);
    rhs(m) = rhs(m) + r*(g1n + g1n_next);

    % solve linear system:
    uint = A1\rhs;

    % augment with boundary values:
    u = [g0n_next; uint; g1n_next];

    % plot results at desired times:
    if mod(n,nplot)==0 || n==nsteps
      ufine = utrue(xfine,t_next);
      plot(x,u,'b.-', 'DisplayName', 'U', xfine, ufine, 'r', 'DisplayName', 'u_{true}')
      title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
      t_next,n,m+2))
      axis([ax bx 0 1])
      error = max(abs(u-utrue(x,t_next)));
      % disp(sprintf('at time t = %9.5e  max error =  %9.5e',t_next,error))
      legend;
      pause(0.01)
    end

    tn = t_next;   % for next time step
  end
  end
