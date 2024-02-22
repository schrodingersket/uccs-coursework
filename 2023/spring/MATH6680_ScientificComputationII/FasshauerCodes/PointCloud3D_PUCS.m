% PointCloud3D_PUCS
% Script that fits a surface to 3D point cloud using partition of
% unity approximation with sparse matrices
% Calls on: CSEvalMatrix
% Uses:     k-D tree package by Guy Shechter
%           from MATLAB Central File Exchange
  % Weight function for global Shepard partition of unity weighting
  wf = @(e,r) r.^4.*(5*spones(r)-4*r);
  % The RBF basis function for local RBF interpolation
  rbf = @(e,r) r.^4.*(5*spones(r)-4*r);
  ep = 1; % Parameter for basis function
  % Parameter for npu-by-npu-by-npu grid of PU cells
  npu = 8;
  % Parameter for npu-by-npu-by-npu grid of PU cells
  neval = 25;
  % Load data points and compute bounding box
  load('Data3D_Bunny3');  N = size(dsites,1);
  bmin = min(dsites,[],1);  bmax = max(dsites,[],1);
  bdim = max(bmax-bmin);
  wep = npu/bdim;
  % Add auxiliary points along normals "inside" and "outside"
  % Find points with nonzero normal vectors and count them
  withnormals = find(normals(:,1)|normals(:,2)|normals(:,3));
  addpoints = length(withnormals);
  % Distance along normal at which to place new points
  delta = bdim/100;
  % Create new points
  dsites(N+1:N+addpoints,:) = ...
       dsites(withnormals,:) + delta*normals(withnormals,:);
  dsites(N+addpoints+1:N+2*addpoints,:) = ...
       dsites(withnormals,:) - delta*normals(withnormals,:);
  % Interpolant is implicit surface, i.e.,
  % "original" points have rhs=0, "inside" rhs=-1, "outside" rhs=1
  rhs = [zeros(N,1); ones(addpoints,1); -ones(addpoints,1)];
  % Compute new bounding box
  bmin = min(dsites,[],1);  bmax = max(dsites,[],1);
  ctrs = dsites;
  % Create neval-by-neval-by-neval equally spaced evaluation
  % locations in bounding box
  xgrid = linspace(bmin(1),bmax(1),neval);
  ygrid = linspace(bmin(2),bmax(2),neval);
  zgrid = linspace(bmin(3),bmax(3),neval);
  [xe,ye,ze] = meshgrid(xgrid,ygrid,zgrid);
  epoints = [xe(:) ye(:) ze(:)];
  % Create npu-by-npu-by-npu equally spaced centers of PU cells
  % in bounding box
  puxgrid = linspace(bmin(1),bmax(1),npu);
  puygrid = linspace(bmin(2),bmax(2),npu);
  puzgrid = linspace(bmin(3),bmax(3),npu);
  [xpu,ypu,zpu] = meshgrid(puxgrid,puygrid,puzgrid);
  cellctrs = [xpu(:) ypu(:) zpu(:)];
  cellradius = 1/wep;
  % Compute Shepard evaluation matrix
  DM_eval = DistanceMatrixCSRBF(epoints,cellctrs,wep);
  SEM = wf(wep,DM_eval);
  SEM = spdiags(1./(SEM*ones(npu^3,1)),0,neval^3,neval^3)*SEM;
  % Build k-D trees for data sites and evaluation points
  [tmp,tmp,datatree] = kdtree(dsites,[]);
  [tmp,tmp,evaltree] = kdtree(epoints,[]);
  Pf = zeros(neval^3,1); % initialize
  for j=1:npu^3
  % Find data sites in cell j
     [pts,dist,idx] = kdrangequery(datatree,...
                     cellctrs(j,:),cellradius);
     if (length(idx) > 0)
        % Build local interpolation matrix for cell j
        DM_data = DistanceMatrixCSRBF(dsites(idx,:),...
                  ctrs(idx,:),ep);
        IM = rbf(ep,DM_data);
        % Find evaluation points in cell j
        [epts,edist,eidx] = kdrangequery(evaltree,...
                            cellctrs(j,:),cellradius);
        % Compute local evaluation matrix
        DM_eval = DistanceMatrixCSRBF(epoints(eidx,:),...
                  ctrs(idx,:),ep);
        EM = rbf(ep,DM_eval);
        % Compute local RBF interpolant
        localfit = EM * (IM\rhs(idx));
        % Accumulate global fit
        Pf(eidx) = Pf(eidx) + localfit.*SEM(eidx,j);
     end
  end
  % Plot data sites with interpolant (zero contour of 3D-fit Pf)
  figure; hold on
  plot3(dsites(1:N,1),dsites(1:N,2),dsites(1:N,3),'bo');
  pfit = patch(isosurface(xe,ye,ze,...
                        reshape(Pf,neval,neval,neval),0));
  isonormals(xe,ye,ze,reshape(Pf,neval,neval,neval),pfit)
  set(pfit,'FaceLighting','gouraud','FaceColor',...
      'red','EdgeColor','none');
  light('Position',[0 0 1],'Style','infinite');
  daspect([1 1 1]); view([0,90]);
  axis([bmin(1) bmax(1) bmin(2) bmax(2) bmin(3) bmax(3)])
  axis off
  hold off

