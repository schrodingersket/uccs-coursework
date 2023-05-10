% ceps = CostEpsilonLRBF(ep,r,rbf,Lrbf)
% Provides the "cost of epsilon" function for LOOCV optimization
% of shape parameter
% Input: ep, values of shape parameter
%        r, dx, Distance and Difference matrices
%        rbf, Lrbf, definition of rbf and its derivative
  function ceps = CostEpsilonLRBF(ep,r,rbf,Lrbf)
  N = size(r,2);
  A = rbf(ep,r);
  rhs = Lrbf(ep,r)';
  invA = pinv(A);
  EF = (invA*rhs)./repmat(diag(invA),1,N);
  ceps = norm(EF(:));
