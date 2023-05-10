% RBFHermite_2D
% Script that performs first-order 2D RBF Hermite interpolation
% Calls on: DistanceMatrix, DifferenceMatrix
  % Define RBF and its derivatives
  rbf = @(e,r) sqrt(1+(e*r).^2);   % MQ RBF
  dxrbf = @(e,r,dx) dx*e^2./sqrt(1+(e*r).^2);
  dyrbf = @(e,r,dy) dy*e^2./sqrt(1+(e*r).^2);
  dxxrbf = @(e,r,dx) e^2*(1+(e*r).^2-(e*dx).^2)./...
                         (1+(e*r).^2).^(3/2);
  dxyrbf = @(e,r,dx,dy) -e^4*dx.*dy./(1+(e*r).^2).^(3/2);
  dyyrbf = @(e,r,dy) e^2*(1+(e*r).^2-(e*dy).^2)./...
                         (1+(e*r).^2).^(3/2);
  ep = 6;
  % Define test function and its derivatives
  tf = @(x,y) (tanh(9*(y-x))+1)/(tanh(9)+1);
  tfDx = @(x,y) 9*(tanh(9*(y-x)).^2-1)/(tanh(9)+1);
  tfDy = @(x,y) 9*(1-tanh(9*(y-x)).^2)/(tanh(9)+1);
  N = 289; gridtype = 'u';
  neval = 40;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute the distance and difference matrices for
  % evaluation matrix
  DM_eval = DistanceMatrix(epoints,ctrs);
  dx_eval = DifferenceMatrix(epoints(:,1),ctrs(:,1));
  dy_eval = Differencematrix(epoints(:,2),ctrs(:,2));
  % Compute the distance and difference matrices for
  % interpolation matrix
  DM_data = DistanceMatrix(dsites,ctrs);
  dx_data = DifferenceMatrix(dsites(:,1),ctrs(:,1));
  dy_data = DifferenceMatrix(dsites(:,2),ctrs(:,2));
  rhs = [tf(dsites(:,1),dsites(:,2)); ...
         tfDx(dsites(:,1),dsites(:,2)); ...
         tfDy(dsites(:,1),dsites(:,2))];
  exact = tf(epoints(:,1),epoints(:,2));
  % Compute blocks for interpolation matrix
  IM = rbf(ep,DM_data);
  DxIM = dxrbf(ep,DM_data,dx_data);
  DyIM = dyrbf(ep,DM_data,dy_data);
  DxxIM = dxxrbf(ep,DM_data,dx_data);
  DxyIM = dxyrbf(ep,DM_data,dx_data,dy_data);
  DyyIM = dyyrbf(ep,DM_data,dy_data);
  % Assemble symmetric interpolation matrix
  IM = [IM -DxIM -DyIM;
        DxIM -DxxIM -DxyIM; 
        DyIM -DxyIM -DyyIM];
  % Compute blocks for evaluation matrix
  EM = rbf(ep,DM_eval);
  DxEM = dxrbf(ep,DM_eval,dx_eval);
  DyEM = dyrbf(ep,DM_eval,dy_eval);
  % Assemble evaluation matrix
  EM = [EM -DxEM -DyEM];
  % RBF Hermite interpolant
  Pf = EM * (IM\rhs);
  % Compute errors on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  fview = [-30,30];  % viewing angles for plot
  caption = 'Hermite interpolant false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  caption = 'Maximum error for Hermite interpolant.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
