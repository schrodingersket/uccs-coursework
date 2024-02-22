% ApproxMLSApprox1D
% Script that performs 1D approximate MLS approximation
% Calls on: DistanceMatrix
  % Laguerre-Gaussians for 1D
  rbf{1} = @(e,r) exp(-(e*r).^2)/sqrt(pi);
  rbf{2} = @(e,r) exp(-(e*r).^2)/sqrt(pi).*(1.5-(e*r).^2);
  rbf{3} = @(e,r) exp(-(e*r).^2)/sqrt(pi).*...
           (1.875-2.5*(e*r).^2+0.5*(e*r).^4);
  D = [2, 4, 6];   % Scale parameters for generating functions
  % Define Franke-like function as testfunction
  f1 = @(x) 0.75*exp(-(9*x-2).^2/4);
  f2 = @(x) 0.75*exp(-(9*x+1).^2/49);
  f3 = @(x) 0.5*exp(-(9*x-7).^2/4);
  f4 = @(x) 0.2*exp(-(9*x-4).^2);
  moll = @(x) 15*exp(-1./(1-4*(x-0.5).^2));
  testfunction = @(x) moll(x).*(f1(x)+f2(x)+f3(x)-f4(x));
  maxlevel = 14;   % number of iterations
  M = 200;   % to create M evaluation points in unit interval
  xe = linspace(0,1,M); epoints = xe(:);
  exact = testfunction(epoints);
  figure; hold on; cword = cellstr(['r- ';'g--';'b: ']); 
  for i=1:length(D)
     for k=1:maxlevel
        N(k) = (2^k+1); ep = (N(k)-1)/sqrt(D(i));
        name = sprintf('Data1D_%du', N(k));  load(name);
        ctrs = dsites;
        % Create vector of function values
        f = testfunction(dsites);
        % Compute evaluation matrix
        DM = DistanceMatrix(epoints,ctrs);
        EM = rbf{i}(ep,DM);
        % Compute approximate MLS approximation
        Pf = EM*f/sqrt(D(i));
        % Compute RMS error on evaluation grid
        rms_err(k) = norm(Pf-exact)/sqrt(M);
     end
     plot(N,rms_err,cword{i});
  end
  set(gca,'XScale','log','YScale','log','Fontsize',14)
  legend('d=0, D=2.0','d=1, D=4.0','d=2, D=6.0',3);
  xlabel('N');  ylabel('Error'); hold off
