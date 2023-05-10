% ceps = CostEpsilon(ep,rbf,r,rhs)
% Implements cost function for optimization of shape parameter 
% via Rippa's LOOCV algorithm
% Example of usage in LOOCV2Dmin.m
  function ceps = CostEpsilon(ep,r,rbf,rhs)
  A = rbf(ep,r);
  invA = pinv(A);
  EF = (invA*rhs)./diag(invA);
  ceps = norm(EF(:),inf);
