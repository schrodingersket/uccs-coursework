% RBFCardinalFunction
% Computes and plots cardinal function for 2D RBF interpolation
% Calls on: DistanceMatrix
  rbf = @(e,r) exp(-(e*r).^2); ep = 5;
  N = 81; gridtype = 'u';
  neval = 80; M = neval^2;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;   % centers coincide with data sites
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  % Compute distance matrix between the data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  % Compute interpolation matrix
  IM = rbf(ep,DM_data);
  % Compute evaluation matrix
  EM = rbf(ep,DM_eval);
  % Compute cardinal functions at evaluation points
  invIM = pinv(IM);
  % centered at datasite(50)
  for j=1:M
     cardvec = (invIM*EM(j,:)')';
     cardfun(j) = cardvec(50);
  end
  figure
  RBFplot = surf(xe,ye,reshape(cardfun,neval,neval));
  set(RBFplot,'FaceColor','interp','EdgeColor','none')
  colormap autumn; view([145 45]); camlight; lighting gouraud
