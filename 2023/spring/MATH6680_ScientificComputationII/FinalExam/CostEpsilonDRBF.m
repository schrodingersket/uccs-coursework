% ceps = CostEpsilonDRBF(ep,r,dx,rbf,dxrbf)
% Provides the "cost of epsilon" function for LOOCV optimization
% of shape parameter
% Input: ep, values of shape parameter
%        r, dx, Distance and Difference matrices
%        rbf, dxrbf, definition of rbf and its derivative
  function ceps = CostEpsilonDRBF(ep,r,dx,rbf,dxrbf)
  N = size(r,2);
  A = rbf(ep,r);   % = A^T since A is symmetric
  rhs = dxrbf(ep,r,dx)';   % A_x^T
  invA = pinv(A);
  EF = (invA*rhs)./repmat(diag(invA),1,N);
  ceps = norm(EF(:));
