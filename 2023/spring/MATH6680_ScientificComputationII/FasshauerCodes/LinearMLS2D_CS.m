% LinearMLS2D_CS
% Script that performs MLS approximation with linear reproduction
% using sparse matrices
% Calls on: LinearScaling2D_CS
  rbf = @(e,r) max(spones(r)-e*r,0).^4.*(4*e*r+spones(r));
  ep = 30; neval = 60;
  % Load data points and rhs
  load('Data2D_DubuqueNE');
  ctrs = dsites;
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute evaluation matrix
  EM = LinearScaling2D_CS(epoints,ctrs,rbf,ep);
  % Compute MLS approximation   (rhs read from data file)
  Pf = EM*rhs;
  figure   % Plot MLS fit
  surfplot = surf(xe,ye,reshape(Pf,neval,neval));
  set(surfplot,'FaceColor','interp','EdgeColor','none')
  view([15,35]); camlight; lighting gouraud; colormap summer
