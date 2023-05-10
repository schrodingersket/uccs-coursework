% ceps = CostEpsilonD2RBF(ep,r,rbf,d2rbf)
% Provides the "cost of epsilon" function for LOOCV optimization
% of shape parameter
% Input: ep, values of shape parameter
%        r, dx, Distance and Difference matrices
%        rbf, dxrbf, definition of rbf and its derivative
  function ceps = CostEpsilonD2RBF(ep,r,rbf,d2rbf)
  N = size(r,2);
  A = rbf(ep,r);
  rhs = d2rbf(ep,r)';
  invA = pinv(A);
  EF = (invA*rhs)./repmat(diag(invA),1,N);
  ceps = norm(EF(:));
