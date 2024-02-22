% [D,x] = DRBF(N,rbf,dxrbf)
% Computes the differentiation matrix D for 1-D derivative
% using Chebyshev points and LOOCV for optimal shape parameter
% Input: N, create N+1 collocation points
%        rbf, dxrbf function handles for rbf and its derivative
% Calls on: DistanceMatrix, DifferenceMatrix
% Requires: CostEpsilonDRBF
  function [D,x] = DRBF(N,rbf,dxrbf)
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)';   % Chebyshev points
  mine = .1; maxe = 10;   % Shape parameter interval
  r = DistanceMatrix(x,x);
  dx = DifferenceMatrix(x,x);
  ep = fminbnd(@(ep) CostEpsilonDRBF(ep,r,dx,rbf,dxrbf),mine,maxe);
  fprintf('Using epsilon = %f\n', ep)
  A = rbf(ep,r);
  Ax = dxrbf(ep,r,dx);
  D = Ax/A;
