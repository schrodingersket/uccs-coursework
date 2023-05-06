figure(1), clf
t = (0:0.1:2*pi)';
noiselevel = 0.1;

% Parametric Fernandez-Guasti squircle
%
r = 1;
s = 0.95;
x = r/(2*s) * sqrt(2 + 2*s*sqrt(2)*cos(t) + s^2*cos(2*t)) - r/(2*s)*sqrt(2 - 2*s*sqrt(2)*cos(t) + s^2*cos(2*t)) + noiselevel * rand(size(t));
y = r/(2*s) * sqrt(2 + 2*s*sqrt(2)*sin(t) - s^2*cos(2*t)) - r/(2*s)*sqrt(2 - 2*s*sqrt(2)*sin(t) - s^2*cos(2*t));

xzero = [x y];

N = length(t); %total number of points in the point cloud

plot(x, y, 'k.')
axis equal

xzero = [xzero(end, :); xzero; xzero(1, :)];

normal = zeros(N, 2);

for j=2:N+1
    data = xzero(j-1:j+1,:);
    normal(j-1,:) = bestfitnormal(data);
end

xzero([1, N + 2], :)=[];
alpha=0.4;
%exterior points

xplus=xzero + alpha*normal;
xminus=xzero - alpha*normal;
plot(xzero(:, 1), xzero(:,2), 'ro', 'DisplayName', 'x_j')
hold on
plot(xplus(:,1), xplus(:, 2), 'b*', 'DisplayName', 'x_j^+')
plot(xminus(:,1), xminus(:, 2), 'gd', 'DisplayName', 'x_j^-')
legend()

axis([-2 2 -2 2])
print('problem_1i.png', '-dpng')

figure(2), clf
plot3(xzero(:, 1), xzero(:, 2), zeros(N, 1), 'ro', 'DisplayName', 'x_j')
hold on
plot3(xplus(:, 1), xplus(:, 2), alpha*ones(N, 1), 'b*', 'DisplayName', 'x_j^+')
plot3(xminus(:, 1), xminus(:, 2), -alpha*ones(N, 1), 'gd', 'DisplayName', 'x_j^-')
axis([-2 2 -2 2 -3*alpha 3*alpha])
legend()
grid on
hold on

eps = 10;
% phi = @(r) exp(-(eps*r).^2); %GA
% phi = @(r) 1./sqrt(1 + (eps*r).^2); %IMQ
phi = @(r) sqrt(1 + (eps*r).^2); %MQ

M = 3*N;
RHS = [zeros(N,1); alpha*ones(N,1); -alpha*ones(N,1)];

A = zeros(M);

Xbig = [xzero; xplus; xminus];

for j = 1:M
    for k = 1:M
        A(j,k) = phi(norm(Xbig(j,:)-Xbig(k,:)));
    end
end

C = A\RHS;

[Xplot, Yplot] = meshgrid(-4:0.05:4, -4:0.05:4);

Zplot = zeros(size(Xplot));
Dist = zeros(M, 1);

rbf = @(Xin) C' * phi(sqrt(sum((Xbig - repmat(Xin, size(Xbig, 1), 1)).^2, 2)));

for j = 1:size(Xplot, 1)
    for k = 1:size(Yplot, 1)
        Zplot(j,k)=rbf([Xplot(j,k) Yplot(j,k)]);
    end
end

surf(Xplot, Yplot, Zplot, 'DisplayName', 'RBF Interpolant')
print('problem_1ii.png', '-dpng')

figure(3), clf
plot(xzero(:, 1), xzero(:,2), 'ro', 'DisplayName', 'x_j')
hold on
plot(xplus(:,1), xplus(:,2), 'b*', 'DisplayName', 'x_j^+')
plot(xminus(:,1), xminus(:,2), 'gd', 'DisplayName', 'x_j^-')
axis([-2 2 -2 2])
v = [-alpha, 0, alpha];
v = [0, 0.001];
contour(Xplot, Yplot, Zplot, v, 'DisplayName', 'RBF Interpolant')
legend()
print('problem_1iii.png', '-dpng')
