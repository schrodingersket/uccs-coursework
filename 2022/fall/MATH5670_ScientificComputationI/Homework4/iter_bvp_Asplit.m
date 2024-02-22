
% Solve 1d BVP obtained by discretizing u''(x) = f(x) 
% with Dirichlet boundary condtions  u(0) = alpha, u(1) = beta
% Discretized using centered difference.
%
% Compare matrix splitting methods:  Jacobi, Gauss-Seidel, and SOR
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)
function [error, rhoG] = iter_bvp_Asplit(method, omega, maxiter, m, ax, bx, alpha, beta, f)
   nplot = 10;          % iterate u is plotted every nplot iterations
    
   h = (bx-ax) / (m+1);
   x = linspace(ax, bx, m+2)';
   
   % determine exact solution to linear system by setting up
   % system Au = f and solving with backslash:
   e = ones(m+2,1);
   
   A = 1/h^2 * spdiags([e -2*e e], [-1 0 1], m+2, m+2);
   A(1,1:2) = [1 0];
   A(m+2,m+1:m+2) = [0 1];
   rhs = f(x);
   rhs(1) = alpha;
   rhs(m+2) = beta;
   ustar = A\rhs;
   
   
   % Decompose A = DA - LA - UA:
   DA = diag(diag(A));   % diagonal part of A
   LA = DA - tril(A);    % strictly lower triangular part of A (negated)
   UA = DA - triu(A);    % strictly upper triangular part of A (negated)
   
   % set up iteration matrix:
   switch method
      case 'Jacobi'
         M = DA;
         N = LA + UA;
      case 'GS'
         M = DA - LA;
         N = UA;
      case 'backward_GS'
         M = DA - UA;
         N = LA;
      case 'SOR'
         % use optimal omega for this problem:
         M = 1/omega * (DA - omega*LA);
         N = 1/omega * ((1-omega)*DA + omega*UA);
   end
   
   
   u = zeros(size(x));   % initial guess for iteration
   
   % figure(1)
   % clf
   % plot(x,ustar,'r')
   % hold on
   % plot(x,u)
   
   error = nan(maxiter+1,1);
   error(1) = max(abs(u-ustar));
   
   
   % Iteration:
   % ----------
   
   
   for iter=1:maxiter

      if strcmp(method, 'naive_SOR')
         % Naive SOR implementation from Exercise 4.1(d):
         uGS = (DA - LA) \ (UA * u + rhs);
         u = u + omega * (uGS - u);
      else
         u = M \ (N*u + rhs);
      end
         
      error(iter+1) = max(abs(u-ustar));
      % if mod(iter,nplot)==0
      %    % plot u every nplot iterations
      %    plot(x,u)
      %    title(sprintf('%s: Iteration %4i, error = %9.3e',...
      %           method,iter,error(iter+1)),'FontSize',15)
      %    drawnow
      %    pause(.1)
      % end
   end
   
   % compute spectral radius of iteration matrix G:
   if strcmp(method, 'naive_SOR')
      rhoG = nan;
   else
      G = M\N;
      rhoG = max(abs(eig(full(G))));
   end
end