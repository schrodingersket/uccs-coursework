warning('off', 'all');

a = [0, 0];   % left boundary
b = [1, 1];   % right boundary

% List of intervals to use for computation in [x, y] form. For instance,
% [4, 5] means that we create 4 intervals in x and 5 intervals in y.
%
intervals = [[9, 9]; [19, 19]; [39, 39]; [59, 59]; [99, 99]];
ntest = length(intervals);

step_sizes_x = zeros(ntest, 1);
step_sizes_y = zeros(ntest, 1);
errors = zeros(ntest, 1);

% f(x,y) function
%
f = @(x, y) 1.25*exp(x + y/2);
laplace_f = @(x, y) 1.25*exp(x + y/2) + 1.25*exp(x + y/2)*(1/2)*(1/2);

% True solution for test problem
%
utrue = @(x, y) exp(x + y/2);

disp(' ')
fprintf('Computed solution on domain [%0.2f, %0.2f] x [%0.2f, %0.2f]:', a(1), b(1), a(2), b(2))
disp(' ')
disp(' ')
fprintf('Grid points      Grid size (hx, hy)      Relative error\n')
fprintf('-------------------------------------------------------\n')
for i = 1:ntest
    [step_sizes_x(i), step_sizes_y(i), errors(i)] = poisson_9pt(b(1), a(1), b(2), a(2), intervals(i, 1), intervals(i, 2), f, laplace_f, utrue, i == ntest);
    fprintf('(%3d, %3d)       (%7.6f, %7.6f)   %10.3e \n', intervals(i, 1) + 1, intervals(i, 2) + 1, step_sizes_x(i), step_sizes_y(i), errors(i));
end

% Estimate order of accuracy from least squares fit:
%
% (see ../fdmbook/matlab/error_loglog.m)
%
Ap = ones(ntest,2);
Ap(:,2) = log(step_sizes_x);
bp = log(errors);
Kp = Ap\bp;
K = Kp(1);
p = Kp(2);
disp(' ')
fprintf('Least squares fit gives E(h) = %g * h^%g',exp(K),p);
disp(' ')
