% PU2D_CS
% Script that performs partition of unity approximation using
% sparse matrices
% Calls on: testfunction2D, DistanceMatrixCSRBF
% Uses:     k-D tree package by Guy Shechter
%           from MATLAB Central File Exchange
  % Weight function for global Shepard partition of unity weighting
  wf = @(e,r) r.^4.*(5*spones(r)-4*r);
  % RBF basis function for local RBF interpolation
  rbf = @(e,r) r.^4.*(5*spones(r)-4*r);
  ep = 0.1; % Parameter for local basis functions
  % Define Franke's function as testfunction
  f1 = @(x,y) 0.75*exp(-((9*x-2).^2+(9*y-2).^2)/4);
  f2 = @(x,y) 0.75*exp(-((9*x+1).^2/49+(9*y+1).^2/10));
  f3 = @(x,y) 0.5*exp(-((9*x-7).^2+(9*y-3).^2)/4);
  f4 = @(x,y) 0.2*exp(-((9*x-4).^2+(9*y-7).^2));
  testfunction = @(x,y) f1(x,y)+f2(x,y)+f3(x,y)-f4(x,y);
  N = 1089; gridtype = 'h';
  % Parameter for npu-by-npu grid of PU cells in unit square
  npu = 16;
  % Parameter for neval-by-neval evaluation grid in unit square
  neval = 40;
  % Load data points
  name = sprintf('Data2D_%d%s',N,gridtype); load(name)
  ctrs = dsites;
  rhs = testfunction(dsites(:,1),dsites(:,2));
  wep = npu; % Parameter for weight function
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Create npu-by-npu equally spaced centers of PU cells in the
  % unit square
  pugrid = linspace(0,1,npu); [xpu,ypu] = meshgrid(pugrid);
  cellctrs = [xpu(:) ypu(:)];
  cellradius = 1/wep;
  % Compute Shepard evaluation matrix
  DM_eval = DistanceMatrixCSRBF(epoints,cellctrs,wep);
  SEM = wf(ep,DM_eval);
  SEM = spdiags(1./(SEM*ones(npu^2,1)),0,neval^2,neval^2)*SEM;
  % Build k-D trees for data sites and evaluation points
  [tmp,tmp,datatree] = kdtree(dsites,[]);
  [tmp,tmp,evaltree] = kdtree(epoints,[]);
  Pf = zeros(neval^2,1);  % initialize
  for j=1:npu^2
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
  % Compute exact solution
  exact = testfunction(epoints(:,1),epoints(:,2));
  % Compute maximum error on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  % Plot interpolant
  fview = [160,20];
  caption = 'PU fit false colored by maximum error.';
  PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
  % Plot maximum error
  caption = 'Maximum error for PU fit.';
  PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)

