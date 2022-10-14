steps = [5, 20, 100, 500];

x0 = 0;
x1 = 1;



disp(' ')
disp('   h     lambda_1   lambda_2   lambda_3   lambda_4   lambda_5')
for i=1:length(steps)
    m = steps(i);
    
    h = (x1 - x0) / m;
    
    e = ones(m, 1);
    A = spdiags([e -2*e e], -1:1, m, m);
    
    alpha = flip(eig(A));
    
    lambda = alpha ./ (h^2);

    disp(sprintf('%1.4f   %5.4f    %5.4f   %5.4f   %5.4f   %5.4f    %5.4f', h, lambda(1), lambda(2), lambda(3), lambda(4), lambda(5)))
end

disp(sprintf('%1.4f   %5.4f    %5.4f   %5.4f   %5.4f   %5.4f    %5.4f', 0, -(1*pi)^2, -(2*pi)^2, -(3*pi)^2, -(4*pi)^2, -(5*pi)^2))