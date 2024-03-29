function [h,k,error, w_avg] = advection_lf_pbc(m, alpha, tfinal, eta)
  %
  % Solve u_t + au_x = 0  on [ax,bx] with periodic boundary conditions,
  % using the Lax-Wendroff method with m interior points.
  %
  % Returns k, h, and the max-norm of the error.
  % This routine can be embedded in a loop on m to test the accuracy,
  % perhaps with calls to error_table and/or error_loglog.
  %
  % From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

  global a
  a = 2;           % advection velocity

  clf              % clear graphics

  ax = 0;
  bx = 1;

  l_inf_0 = 0;
  l_inf_tfinal = 0;

  h = (bx-ax)/(m+1);         % h = delta x
  k = alpha*h;                  % time step
  nu = a*k/h;                % Courant number
  x = linspace(ax,bx,m+2)';  % note x(1)=0 and x(m+2)=1
  % With periodic BC's there are m+1 unknowns u(2:m+2)
  I = 2:(m+2);   % indices of unknowns

  nsteps = round(tfinal / k);    % number of time steps
  nplot = 20;       % plot solution every nplot time steps
  % (set nplot=2 to plot every 2 time steps, etc.)
  nplot = nsteps;  % only plot at final time

  if abs(k*nsteps - tfinal) > 1e-5
    % The last step won't go exactly to tfinal.
    disp(' ')
    disp(sprintf('WARNING *** k does not divide tfinal, k = %9.5e',k))
    disp(' ')
  end

  % initial conditions:
  tinit = 0;
  tn = tinit;
  u0 = eta(x);
  u = u0;

  % periodic boundary conditions:
  u(1) = u(m+2);   % copy value from rightmost unknown to ghost cell on left
  u(m+3) = u(2);   % copy value from leftmost unknown to ghost cell on right

  % initial data on fine grid for plotting:
  xfine = linspace(ax,bx,1001);
  ufine = utrue(xfine,0, eta);

  % plot initial data:
  plot(x,u0,'b.-', xfine,ufine,'r')
  axis([0 1 -1 1.2])
  legend('computed','true')
  title('Initial data at time = 0')

%   input('Hit <return> to continue  ');

  % main time-stepping loop:
  u_prev = utrue(x, tn - k, eta)(I);
  [max_val, l_inf_0] = max(u);
  disp(sprintf('Initial wave max is at x=%0.2f', ax + l_inf_0 * h));
  for n = 1:nsteps
    tnp = tn + k;   % = t_{n+1}

    % leapfrog:
    u_current = u(I); % store current solution to save u_prev
    u(I) = u_prev - nu*(u(I+1) - u(I-1));
    u_prev = u_current; % set previous timestep to the solution before the current step

    % periodic boundary conditions:
    u(1) = u(m+2);   % copy value from rightmost unknown to ghost cell on left
    u(m+3) = u(2);   % copy value from leftmost unknown to ghost cell on right

    % plot results at desired times:
    if mod(n,nplot)==0 || n==nsteps
      uint = u(1:m+2);  % points on the interval (drop ghost cell on right)
      ufine = utrue(xfine,tnp, eta);
      plot(x,uint,'b.-', xfine,ufine,'r')
      axis([0 1 -1 1.2])
      title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
      tnp,n,m+1))
      error = max(abs(uint-utrue(x,tnp, eta)));
      disp(sprintf('at time t = %9.5e  max error =  %9.5e',tnp,error))
      if n < nsteps, input('Hit <return> to continue  '); end;
    end

    tn = tnp;   % for next time step
  end
  [max_val, l_inf_tfinal] = max(u);
  disp(sprintf('Final wave max is at x=%0.2f', ax + l_inf_tfinal * h));

  w_avg = h*(l_inf_tfinal - l_inf_0) / (tfinal - tinit);
  disp(sprintf('Estimated group velocity is %0.2f', w_avg));

  %--------------------------------------------------------

  function utrue = utrue(x,t, eta)
    % true solution for comparison
    global a

    % For periodic BC's, we need the periodic extension of eta(x).
    % Map x-a*t back to unit interval:

    xat = rem(x - a*t, 1);
    ineg = find(xat<0);
    xat(ineg) = xat(ineg) + 1;
    utrue = eta(xat);
    return