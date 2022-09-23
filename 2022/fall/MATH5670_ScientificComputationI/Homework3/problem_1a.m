a = 0;   % left boundary
b = 1;   % right boundary

intervals = [19, 49, 99, 199, 499];
ntest = length(intervals);
step_size = zeros(ntest, 1);
errors = zeros(ntest, 1);

% f(x,y) function
%  
f = @(x, y) 1.25*exp(x + y/2);

% True solution for test problem
%
utrue = @(x, y) exp(x + y/2);

fprintf('Grid points      Grid size (h)           Relative error\n')
fprintf('-------------------------------------------------------\n')
for i = 1:length(intervals)
    [step_size(i), errors(i)] = poisson(a, b, intervals(i), f, utrue);

    fprintf('%3d            %10.6f               %10.3e \n', intervals(i) + 1, step_size(i), errors(i));
end

% Estimate order of accuracy from least squares fit:
%
% (see ../fdmbook/matlab/error_loglog.m)
%
Ap = ones(ntest,2);
Ap(:,2) = log(step_size);
bp = log(errors);
Kp = Ap\bp;
K = Kp(1);
p = Kp(2);
disp(' ')
disp(sprintf('Least squares fit gives E(h) = %g * h^%g',exp(K),p))
disp(' ')