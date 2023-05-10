% [D2,x] = D2RBF(N,rbf,d2rbf)
% Computes the second-order differentiation matrix D2 for 1-D
% derivative using Chebyshev points and LOOCV for optimal epsilon
% Input: N, number of points -1
%        rbf, d2rbf, function handles for rbf and its derivative
% Calls on: DistanceMatrix, DifferenceMatrix
% Requires: CostEpsilonD2RBF
  function [D2,x] = D2RBF(N,rbf,d2rbf)
  if N==0, D2=0; x=1; return, end
  x = cos(pi*(0:N)/N)';   % Chebyshev points
  mine = .1; maxe = 10;   % Shape parameter interval
  r = DistanceMatrix(x,x);
  ep = fminbnd(@(ep) CostEpsilonD2RBF(ep,r,rbf,d2rbf),mine,maxe);
  fprintf('Using epsilon = %f\n', ep)
  A = rbf(ep,r);
  AD2 = d2rbf(ep,r);
  D2 = AD2/A;
