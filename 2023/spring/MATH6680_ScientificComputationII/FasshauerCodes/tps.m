% rbf = tps(e,r)
% Defines thin plate spline RBF
function rbf = tps(e,r) 
rbf = zeros(size(r));
nz = find(r~=0);   % to deal with singularity at origin
rbf(nz) = (e*r(nz)).^2.*log(e*r(nz));
