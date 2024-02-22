% p12.m - accuracy of Chebyshev spectral differentiation
% (compare p7.m)
% Compute second derivatives for various values of N:
Nmax = 50; E = zeros(4,Nmax);
for N = 1:Nmax;
  [D, chebx] = cheb(N); % differentiation matrix
  D2 = D^2;

  v = @(x) abs(x).^3; % 3rd deriv in BV
  [xx, vv, w2] = chebfft2(v(chebx));
  E(1,N) = norm(w2 - D2 * v(chebx),inf);

  v = @(x) exp(-x.^(-2)); % C-infinity
  [xx, vv, w2] = chebfft2(v(chebx));
  E(2,N) = norm(w2 - D2 * v(chebx),inf);

  v = @(x) 1./(1+x.^2); % analytic in [-1,1]
  [xx, vv, w2] = chebfft2(v(chebx));
  E(3,N) = norm(w2 - D2 * v(chebx),inf);

  v = @(x) x.^10; % polynomial
  [xx, vv, w2] = chebfft2(v(chebx));
  E(4,N) = norm(w2 - (D2 * v(chebx)), inf);
end

% Plot results:
titles = {'|x^3|','exp(-x^{-2})','1/(1+x^2)','x^{10}'}; clf
for iplot = 1:size(E)(1)
  subplot(2,2,iplot)
  semilogy(1:Nmax,E(iplot,:),'.','markersize',6)
  line(1:Nmax,E(iplot,:), 'linewidth',.8)
  axis([0 Nmax 1e-16 1e3]), grid on
  set(gca,'xtick',0:10:Nmax,'ytick',(10).^(-15:5:0))
  xlabel N, ylabel error, title(titles(iplot))
end

print('problem_7', '-dpng')
