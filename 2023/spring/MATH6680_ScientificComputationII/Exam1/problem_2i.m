u = @(z) (1 + exp(z/sqrt(6)))^(-2);
u_prime = @(z) -2/sqrt(6) * exp(z/sqrt(6)) * (1 + exp(z/sqrt(6)))^(-3);

L = 20

function dphidt = fisher(t,y)
    c = 5/sqrt(6);
    dphidt = [y(2); -c * y(2) - y(1)*(1 - y(1))];
end

[t, y] = ode45(@fisher, [-L L], [u(-L); u_prime(-L)]);

plot(t,y(:,1),'-o')
title('Fisher Equation ODE Solution (c = 5/sqrt(6))');
xlabel('z');
ylabel('\phi');

print('-dpng', 'problem_2i')