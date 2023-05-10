% ML_HermiteLaplaceCSRBF2D
% Script that performs symmetric multilevel RBF collocation
% using sparse matrices
% Calls on: DistanceMatrixCSRBF
  % Wendland C6 RBF, its Laplacian and double Laplacian
  rbf = @(e,r) r.^8.*(66*spones(r)-154*r+121*r.^2-32*r.^3);
  Lrbf = @(e,r) 44*e^2*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
  L2rbf = @(e,r) 1056*e^4*r.^4.*...
                 (105*spones(r)-483*r+679*r.^2-297*r.^3);
  % Exact solution and its Laplacian for test problem
  u = @(x,y) sin(pi*x).*cos(pi*y/2);
  Lu = @(x,y) -1.25*pi^2*sin(pi*x).*cos(pi*y/2);
  K = 6; neval = 40; gridtype = 'h';
  ep = 0.5*2.^[0:K-1];
  % Create neval-by-neval equally spaced evaluation locations
  % in the unit square
  grid = linspace(0,1,neval); [xe,ye] = meshgrid(grid);
  epoints = [xe(:) ye(:)];
  % Compute exact solution
  exact = u(epoints(:,1),epoints(:,2));
  Rf_old = zeros(17,1);
  for k=1:K
     N1 = (2^k+1)^2;    N2 = (2^(k+1)+1)^2;
     name1 = sprintf('Data2D_%d%s',N1,gridtype);
     name2 = sprintf('Data2D_%d%s',N2,gridtype);
     load(name2)
     % Additional boundary points for residual evaluation
     sn = sqrt(N2); bdylin = linspace(0,1,sn)';
     bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
     bdyres = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1); ...
          flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
     intres = dsites;
     load(name1); intdata = dsites;
     % Additional boundary points
     sn = sqrt(N1); bdylin = linspace(0,1,sn)';
     bdy0 = zeros(sn-1,1); bdy1 = ones(sn-1,1);
     bdydata = [bdylin(1:end-1) bdy0; bdy1 bdylin(1:end-1); ...
          flipud(bdylin(2:end)) bdy1; bdy0 flipud(bdylin(2:end))];
     bdyctrs{k} = bdydata;
     intctrs{k} = intdata;
     % Compute new right-hand side (= residual)
     Tf = [Lu(intdata(:,1),intdata(:,2)); ...
           sin(pi*bdydata(1:sn-1,1)); zeros(3*(sn-1),1)];
     rhs = Tf - Rf_old;
     % Compute blocks for collocation matrix
     DM_IIdata = DistanceMatrixCSRBF(intdata,intctrs{k},ep(k));
     DM_IBdata = DistanceMatrixCSRBF(intdata,bdyctrs{k},ep(k));
     DM_BIdata = DistanceMatrixCSRBF(bdydata,intctrs{k},ep(k));
     DM_BBdata = DistanceMatrixCSRBF(bdydata,bdyctrs{k},ep(k));
     LLCM = L2rbf(ep(k),DM_IIdata);
     LBCM = Lrbf(ep(k),DM_IBdata);
     BLCM = Lrbf(ep(k),DM_BIdata);
     BBCM = rbf(ep(k),DM_BBdata);
     CM = [LLCM LBCM; BLCM BBCM];
     % Compute coefficients for RBF solution of detail level
     coef{k} = CM\rhs;
     if (k < K)
        % based on the distances between the next finer
        % points (respoints) and centers
        for j=1:k
           DM_IIres = DistanceMatrixCSRBF(intres,intctrs{j},ep(j));
           DM_IBres = DistanceMatrixCSRBF(intres,bdyctrs{j},ep(j));
           DM_BIres = DistanceMatrixCSRBF(bdyres,intctrs{j},ep(j));
           DM_BBres = DistanceMatrixCSRBF(bdyres,bdyctrs{j},ep(j));
           LLRM = L2rbf(ep(j),DM_IIres);
           LBRM = Lrbf(ep(j),DM_IBres);
           BLRM = Lrbf(ep(j),DM_BIres);
           BBRM = rbf(ep(j),DM_BBres);
           RM{j} = [LLRM LBRM; BLRM BBRM];
        end
        % Evaluate RBF approximation (sum of all previous fits,
        % but evaluated on current grid)
        Rf = zeros(N2+4*sqrt(N2)-4,1);
        for j=1:k
           Rf = Rf + RM{j}*coef{j};
        end
        Rf_old = Rf;
     end
     % Compute evaluation matrix
     DM_inteval = DistanceMatrixCSRBF(epoints,intctrs{k},ep(k));
     DM_bdyeval = DistanceMatrixCSRBF(epoints,bdyctrs{k},ep(k));
     LEM = Lrbf(ep(k),DM_inteval);
     BEM = rbf(ep(k),DM_bdyeval);
     EM = [LEM BEM];
     % Evaluate RBF approximation
     Pf = EM*coef{k};
     if (k > 1)
        Pf = Pf_old + Pf;
     end
     Pf_old = Pf;
     % Compute maximum error on evaluation grid
     maxerr = norm(Pf-exact,inf);
     rms_err = norm(Pf-exact)/neval;
     fprintf('RMS error:     %e\n', rms_err)
     fprintf('Maximum error: %e\n', maxerr)
     if (k > 1)
        max_rate = log(maxerr_old/maxerr)/log(2);
        rms_rate = log(rms_err_old/rms_err)/log(2);
        fprintf('RMS rate:      %f\n', rms_rate)
        fprintf('Maxerror rate: %f\n', max_rate)
     end
     maxerr_old = maxerr;  rms_err_old = rms_err;
     % Plot collocation solution
     fview = [-30,30];
     caption = ['ML RBF solution for level ',num2str(k),...
         ' false colored by maximum error'];
     PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
     caption = ['Maximum error for ML RBF solution for level ',...
         num2str(k)];
     PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption);
  end
