% DistanceMatrixFit
% Script that uses Euclidean distance matrices to perform
% scattered data interpolation for arbitrary space dimensions
% Calls on: DistanceMatrix, MakeSDGrid, testfunction
% Uses:     haltonseq (written by Daniel Dougherty from Matlab
%                      Central File Exchange)
  s = 2;
  k = 3; N = (2^k+1)^s;
  neval = 10; M = neval^s;
  % Use Halton points as data sites and centers
  dsites = haltonseq(N,s);
  ctrs = dsites;
  % Create neval^s equally spaced evaluation locations in the
  % s-dimensional unit cube
  epoints = MakeSDGrid(s,neval);
  % Create right-hand side vector,
  % i.e., evaluate the test function at the data sites
  rhs = testfunction(s,dsites);
  % Compute distance matrix for the data sites and centers
  IM = DistanceMatrix(dsites,ctrs);
  % Compute distance matrix for evaluation points and centers
  EM = DistanceMatrix(epoints,ctrs);
  % Evaluate the interpolant on evaluation points
  % (evaluation matrix * solution of interpolation system)
  Pf = EM * (IM\rhs);
  % Compute exact solution,
  % i.e., evaluate test function on evaluation points
  exact = testfunction(s,epoints);
  % Compute maximum and RMS errors on evaluation grid
  maxerr = norm(Pf-exact,inf);
  rms_err = norm(Pf-exact)/sqrt(M);
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  switch s
     case 1
        plot(epoints, Pf)
        set(gca,'Fontsize',14)
        xlabel('x','FontSize',14);
        ylabel('y','FontSize',14,'Rotation',0);
        figure; plot(epoints, abs(Pf-exact))
        set(gca,'Fontsize',14)
        xlabel('x','FontSize',14);
        ylabel('Error','FontSize',14);
     case 2
        fview = [-30,30];
        xe = reshape(epoints(:,2),neval,neval);
        ye = reshape(epoints(:,1),neval,neval);
        caption = ['Distance matrix fit ',...
             'false colored by maximum error.'];
        PlotSurf(xe,ye,Pf,neval,exact,maxerr,fview,caption);
        caption = 'Maximum error for distance matrix fit.';
        PlotError2D(xe,ye,Pf,exact,maxerr,neval,fview,caption)
     case 3
        xe = reshape(epoints(:,2),neval,neval,neval);
        ye = reshape(epoints(:,1),neval,neval,neval);
        ze = reshape(epoints(:,3),neval,neval,neval);
        xslice = .25:.25:1; yslice = 1; zslice = [0,0.5];
        caption = 'Slice plot of distance matrix fit.';
        PlotSlices(xe,ye,ze,Pf,neval,xslice,yslice,zslice,caption);
        caption = ['Slice plot of absolute error ',...
            'for distance matrix fit.'];
        PlotErrorSlices(xe,ye,ze,Pf,exact,neval,...
            xslice,yslice,zslice,caption);
     otherwise
        disp('Cannot display plots for s>3')   
  end
