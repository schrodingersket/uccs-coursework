rng(0)

clf

N=10;
x = rand(1,N);
y = rand(1,N);
z = (3-2*x-5*y)/4; % Equation of the plane containing 
                   % (x,y,z) points is 2*x+5*y+4*z=3

% Turn X, Y, Z into column vectors
Xcol = x(:); 
Ycol = y(:); 
Zcol = z(:); 
Const = ones(size(Xcol)); % Vector of the constant term in the RHS
n = [Xcol Ycol Zcol]\Const; % Find the coefficients
XCoeff = n(1); % X coefficient
YCoeff = n(2); % Y coefficient
ZCoeff = n(3); % Z coefficient

% Using the above variables, z = XCoeff * x + YCoeff * y + CCoeff
L=plot3(x,y,z,'ro'); % Plot the original data points
set(L,'Markersize',2*get(L,'Markersize')) % Making the circle markers larger
set(L,'Markerfacecolor','r') % Filling in the markers
hold on
[xx, yy]=meshgrid(min(x):0.1:max(x),min(y):0.1:max(y)); % Generating a regular grid for plotting
zz = (1 - XCoeff * xx - YCoeff * yy)/ZCoeff;
surf(xx,yy,zz) % Plotting the surface
title(sprintf('Plotting plane (%f)*x + (%f)*y + (%f)*z = 1',XCoeff, YCoeff, ZCoeff))

print('problem_1_normal_ex', '-dpng')