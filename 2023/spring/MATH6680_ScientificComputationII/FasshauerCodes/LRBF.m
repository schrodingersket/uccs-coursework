% [L,x,y] = LRBF(N,rbf,Lrbf)
% Computes the Laplacian differentiation matrix L for 2-D
% derivatives using Chebyshev points and LOOCV for optimal epsilon
% Input: N number of points -1
%        rbf, Lrbf, function handles for rbf and its derivative
% Calls on: DistanceMatrix
% Requires: CostEpsilonLRBF
  function [L,x,y] = LRBF(N,rbf,Lrbf)
  if N==0, L=0; x=1; return, end
  x = cos(pi*(0:N)/N)';   % Chebyshev points
  y = x; [xx,yy] = meshgrid(x,y);
  % Stretch 2D grids to 1D vectors and put in one array
  points = [xx(:) yy(:)];
  mine = .1; maxe = 10;   % Shape parameter interval
  r = DistanceMatrix(points,points);
  ep = fminbnd(@(ep) CostEpsilonLRBF(ep,r,rbf,Lrbf),mine,maxe);
  fprintf('Using epsilon = %f\n', ep)
  A = rbf(ep,r);
  AL = Lrbf(ep,r);
  L = AL/A;
