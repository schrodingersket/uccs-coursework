% Phi = LinearScaling2D_CS(epoints,ctrs,rbf,ep)
% Forms a sparse matrix of scaled generating functions for MLS
% approximation with linear reproduction.
% Uses:         k-D tree package by Guy Shechter from
%               MATLAB Central File Exchange
  function Phi = LinearScaling2D_CS(epoints,ctrs,rbf,ep)
  [N,s] = size(epoints); [M,s] = size(ctrs);
  alpha = [0 0; 1 0; 0 1; 1 1; 2 0; 0 2];
  % Build k-D tree for centers
  [tmp,tmp,Tree] = kdtree(ctrs,[]);
  % For each eval. point, find centers whose support overlap it
  support = 1/ep; mu = zeros(6);
  % Modify the following line for optimum performance
  veclength = round(support*N*M/4);
  rowidx = zeros(1,veclength); colidx = zeros(1,veclength);
  validx = zeros(1,veclength);
  istart = 1; iend = 0;
  for i = 1:N
     [pts,dist,idx] = kdrangequery(Tree,epoints(i,:),support);
     newlen = length(idx);
     % Vector of basis functions
     Phi_i = rbf(ep,dist');
     % Compute all 6 moments for i-th evaluation point
     for j=1:6
        x_to_alpha = 1;
        for coord=1:s
           x_to_alpha = x_to_alpha .*(ctrs(idx,coord)-...
             repmat(epoints(i,coord),newlen,1)).^alpha(j,coord);
        end
        mu(j) = Phi_i*x_to_alpha;
     end
     L1=(mu(4)^2-mu(5)*mu(6)); L2=(mu(2)*mu(6)-mu(3)*mu(4));
     L3=(mu(5)*mu(3)-mu(2)*mu(4));
     scaling = L1*repmat(1,newlen,1) + ...
           L2*(ctrs(idx,1)-repmat(epoints(i,1),newlen,1))+...
           L3*(ctrs(idx,2)-repmat(epoints(i,2),newlen,1));
     denom = mu(2)^2*mu(6)+mu(5)*mu(3)^2-mu(1)*mu(5)*mu(6)-...
             2*mu(2)*mu(3)*mu(4)+mu(1)*mu(4)^2;
     if (denom ~= 0)
        scaling = scaling/denom;
        iend = iend + newlen;
        rowidx(istart:iend) = repmat(i,1,newlen);
        colidx(istart:iend) = idx';
        validx(istart:iend) = Phi_i.*scaling';
        istart = istart + newlen;
     end
  end
  filled = find(rowidx); % only those actually filled
  Phi = sparse(rowidx(filled),colidx(filled),validx(filled),N,M);
  % Free memory
  clear rowidx colidx validx; kdtree([],[],Tree);
