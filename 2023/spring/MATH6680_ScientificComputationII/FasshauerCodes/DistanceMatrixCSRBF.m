% DM = DistanceMatrixCSRBF(dsites,ctrs,ep)
% Forms the distance matrix of two sets of points in R^s
% for compactly supported radial basis functions, i.e.,
%      DM(i,j) = || datasite_i - center_j ||_2.
% The CSRBF used with this code must be given in shifted form
% rbf2(u) = rbf(r), u=1-e*r.
% For example, the Wendland C2
% rbf = @(e,r) max(1-e*r,0).^4.*(4*e*r+1);
% becomes
% rbf2 = @(u) u.^4.*(4*u+5);
% Input
%   dsites: Nxs matrix representing a set of N data sites
%              in R^s (i.e., each row contains one
%              s-dimensional point)
%   ctrs:   Mxs matrix representing a set of M centers for
%              RBFs in R^s (also one center per row)
%   ep:        determines size of support of basis function.
%              Small ep yields wide function,
%              i.e., supportsize = 1/ep
% Output
%   DM:     NxM SPARSE matrix that contains the Euclidean
%              u-distance (u=1-e*r) between the i-th data
%              site and the j-th center in the i,j position
% Uses:     k-D tree package by Guy Shechter from
%              MATLAB Central File Exchange
  function DM = DistanceMatrixCSRBF(dsites,ctrs,ep)
  N = size(dsites,1);  M = size(ctrs,1);
  % Build k-D tree for data sites
  % For each center (basis function), find the data sites
  % in its support along with u-distance
  support = 1/ep;
  nzmax = 25*N; rowidx = zeros(1,nzmax); colidx = zeros(1,nzmax);
  validx = zeros(1,nzmax); istart = 1; iend = 0;
  if M > N  % faster if more centers than data sites
     [tmp,tmp,Tree] = kdtree(ctrs,[]);
     for i = 1:N
        [pts,dist,idx] = kdrangequery(Tree,dsites(i,:),support);
        newentries = length(idx);
        iend = iend + newentries;
        rowidx(istart:iend) = repmat(i,1,newentries);
        colidx(istart:iend) = idx';
        validx(istart:iend) = 1-ep*dist';
        istart = istart + newentries;
     end
  else
     [tmp,tmp,Tree] = kdtree(dsites,[]);
     for j = 1:M
        [pts,dist,idx] = kdrangequery(Tree,ctrs(j,:),support);
        newentries = length(idx);
        iend = iend + newentries;
        rowidx(istart:iend) = idx';
        colidx(istart:iend) = repmat(j,1,newentries);
        validx(istart:iend) = 1-ep*dist';
        istart = istart + newentries;
     end
  end
  idx = find(rowidx);
  DM = sparse(rowidx(idx),colidx(idx),validx(idx),N,M);
  % Free the k-D Tree from memory.
  kdtree([],[],Tree);
