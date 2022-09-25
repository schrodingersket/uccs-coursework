order = 4;

% Point of interest
%
x_bar = 1;

% Function definitions
%
u = @(x) sin(2*x);
uxx = @(x) -4*sin(2*x);
ux6 = @(x) -(2^6)*sin(2*x);
dd_4 = @(x, h) (1/h^2) * ((-1/12)*u(x-2*h) + (4/3)*u(x-h) + (-5/2)*u(x) + (4/3)*u(x+h) + (-1/12)*u(x+2*h));

uxx_true = uxx(x_bar);
hvals = logspace(-1, -4, 13);
E_DD4u = [];
E_DD4u_predicted = [];

% table headings:
disp(' ')
disp('       h           Error [DD4u]    Predicted Error [DD4u]')

for i=1:length(hvals)
   h = hvals(i);

   % Fourth-order approximation to u''(x_bar):
   %
   DD4u = dd_4(x_bar, h);

   % Error:
   %
   E_DD4u(i) = DD4u - uxx_true;

   % Predicted error:
   %
   E_DD4u_predicted(i) = (-1/90) * h^4 * ux6(x_bar);

   % Print line of table:
   %
   disp(sprintf('%13.4e   %13.4e   %13.4e', h, E_DD4u(i), E_DD4u_predicted(i)))
end

% Plot errors
%
clf;
axis([5e-4 .2 1e-12 1]);
plt = loglog(hvals, abs(E_DD4u), 'o-');
grid on;
title('Step size vs. Absolute Error')
xlabel('Step size (h)')
ylabel('Error (L1)')
waitfor(plt);