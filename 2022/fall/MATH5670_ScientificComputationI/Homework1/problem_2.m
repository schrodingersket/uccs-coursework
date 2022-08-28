#!/usr/bin/env /usr/bin/octave

order = 4;

% Point of interest
%
x_bar = 1;

% Function definitions
%
u = @(x) sin(2*x);
upp = @(x) -4*sin(2*x);
dd_4 = @(x, h) (1/h^2) * ((-1/12)*u(x-2*h) + (4/3)*u(x-h) + (-5/2)*u(x) + (4/3)*u(x+h) + (-1/12)*u(x+2*h));

upp_true = upp(x_bar);
hvals = logspace(-1, -4, 13);
E_DD4u = [];

% table headings:
disp(' ')
disp('       h              DD4u')

for i=1:length(hvals)
   h = hvals(i);

   % Approximation to u''(1):
   %
   DD4u = dd_4(1, h);

   % Error:
   %
   E_DD4u(i) = DD4u - upp_true;

   % print line of table:
   disp(sprintf('%13.4e   %13.4e', h, E_DD4u(i)))
end

% Plot errors
%
clf;
axis([5e-4 .2 1e-12 1]);
plt = loglog(hvals, abs(E_DD4u), 'o-');
waitfor(plt);