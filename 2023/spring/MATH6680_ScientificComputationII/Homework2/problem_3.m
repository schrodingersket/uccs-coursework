D5 = cheb(5);
D20 = cheb(20);

fprintf('||(D5)^6||_2 = %0.10f\n', norm(D5^6, 2))
fprintf('||(D20)^21||_2 = %0.10f\n', norm(D20^21, 2))