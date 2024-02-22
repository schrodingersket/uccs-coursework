% ML_CSRBF3D
% Script that performs multilevel RBF Interpolation using
% sparse matrices
% Calls on: DistancMatrixCSRBF
  % Wendland C4
  rbf = @(e,r) r.^6.*(35*r.^2-88*r+56*spones(r));
  % Number of levels with epsilons for stationary interpolation
  K = 4; ep = 0.7*2.^[0:K-1];
  testfunction = @(x,y,z) 64*x.*(1-x).*y.*(1-y).*z.*(1-z);
  gridtype = 'u'; % Type of data points: 'u'=uniform
  neval = 10; M = neval^3;
  grid=linspace(0,1,neval); [xe,ye,ze]=meshgrid(grid);
  epoints=[xe(:) ye(:) ze(:)];
  exact = testfunction(epoints(:,1),epoints(:,2),epoints(:,3));
  isomin = 0.1; isomax = 1; isostep = .1;
  xslice = .25:.25:1; yslice = 1; zslice = [0,0.5];
  Rf_old = zeros(27,1);    % initialize
  for k=1:K
     N1 = (2^k+1)^3;    N2 = (2^(k+1)+1)^3;
     name1 = sprintf('Data3D_%d%s',N1,gridtype);
     name2 = sprintf('Data3D_%d%s',N2,gridtype);
     load(name2);     respoints = dsites;
     load(name1);     ctrs{k} = dsites;
     % Compute right-hand side (= residual)
     Tf = testfunction(dsites(:,1),dsites(:,2),dsites(:,3));
     rhs = Tf - Rf_old;
     DM_data = DistanceMatrixCSRBF(dsites,ctrs{k},ep(k));
     IM = rbf(ep(k),DM_data);
     % Compute coefficients for RBF interpolant to detail level
     coef{k} = IM\rhs;
     if (k < K)
        % Compute - for all levels - evaluation matrices for
        % residuals directly
        for j=1:k
           DM_res = DistanceMatrixCSRBF(respoints,ctrs{j},ep(j));
           RM{j} = rbf(ep(j),DM_res);
        end
        % Evaluate RBF interpolant (sum of all previous fits
        % evaluated on current grid)
        Rf = zeros(N2,1);
        for j=1:k
           Rf = Rf + RM{j}*coef{j};
        end
        Rf_old = Rf;
     end
     DM_eval = DistanceMatrixCSRBF(epoints,ctrs{k},ep(k));
     EM = rbf(ep(k),DM_eval);
     Pf = EM*coef{k};
     if (k > 1)
        Pf = Pf_old + Pf;
     end
     Pf_old = Pf;
     maxerr = norm(Pf-exact,inf);
     rms_err = norm(Pf-exact)/sqrt(M);
     fprintf('RMS error:     %e\n', rms_err)
     if (k > 1)
        rms_rate = log(rms_err_old/rms_err)/log(2);
        fprintf('RMS rate:      %f\n', rms_rate)
     end
     rms_err_old = rms_err;
     % Plot data sites
     figure
     plot3(dsites(:,1),dsites(:,2),dsites(:,3),'bo');
     caption = ['Isosurfaces of ML CSRBF interpolant for level '...
          num2str(k)];
     PlotIsosurf(xe,ye,ze,Pf,neval,exact,maxerr,isomin,...
         isostep,isomax,caption);
     caption = ['Slice plot of ML CSRBF interpolant for level '...
          num2str(k)];
     PlotSlices(xe,ye,ze,Pf,neval,xslice,yslice,zslice,caption);
     caption = ['Slice plot of absolute error',... 
         ' for ML RBF interpolant for level ' num2str(k)];
     PlotErrorSlices(xe,ye,ze,Pf,exact,neval,...
         xslice,yslice,zslice,caption);
  end
