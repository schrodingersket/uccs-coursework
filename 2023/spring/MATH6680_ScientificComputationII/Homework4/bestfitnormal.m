function out=bestfitnormal(data)

[N, dim]=size(data);

Const = ones(N,1); % Vector of the constant term in the RHS
out = data\Const; % Find the coefficients

switch dim
    case 2
        x=data(:,1); y=data(:,2); L=plot(x,y,'ro'); % Plot the original data points
        % set(L,'Markersize',2*get(L,'Markersize')) % Making the circle markers larger
        % set(L,'Markerfacecolor','r') % Filling in the markers
        % hold on
        % xx=min(x):0.1:max(x);
        % yy=(1-out(1)*xx)/out(2);
        % line(xx,yy) % plotting the line
        % title(sprintf('Plotting line (%f)*x+(%f)*y=1',out(1), out(2)))
        % hold off
    case 3
        x=data(:,1); y=data(:,2); z=data(:,3);
        L=plot3(x,y,z,'ro'); % Plot the original data points
        % set(L,'Markersize',2*get(L,'Markersize')) % Making the circle markers larger
        % set(L,'Markerfacecolor','r') % Filling in the markers
        % hold on
        % [xx, yy]=meshgrid(min(x):0.1:max(x),min(y):0.1:max(y)); % Generating a regular grid for plotting
        % zz = (1 - out(1) * xx - out(2) * yy)/out(3);
        % surf(xx,yy,zz) % Plotting the surface
        % title(sprintf('Plotting plane (%f)*x+(%f)*y+(%f)*z=1',out(1), out(2), out(3)))
        % hold off
    otherwise
        fprintf('The data is not 2D or 3D')

end
