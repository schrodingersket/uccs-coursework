% CHEBFFT Chebyshev differentiation via FFT. Simple, not optimal.
% If v is complex, delete "real" commands.
function [x, v, w] = chebfft2(v)
  N = length(v)-1; 
  if N==0
    x=0;
    w=0;
  return, end
  x = cos((0:N)'*pi/N);
  ii = 0:N-1;
  v = v(:); V = [v; flipud(v(2:N))]; % transform x -> theta
  U = real(fft(V));
  W1 = real(ifft((1j)^1 .* ([ii 0 1-N:-1].^1)' .* U));
  W2 = real(ifft((1j)^2 .* ([ii 0 1-N:-1].^2)' .* U));
  w = zeros(N+1,1);

  % transform theta -> x

  % Interior points (Trefethen 8.7)
  %
  w(2:N) = W1(2:N) .* (-x(2:N) ./ (1 - x(2:N).^2).^(3/2)) + W2(2:N) ./ (1 - x(2:N).^2);

  % End points (Trefethen 8.7)
  %
  w(1) = (2/3) * sum((ii'.^2 + ii'.^4).*U(ii+1))/N + .5*N*U(N+1);

  w(N+1) = (2/3) * sum((-1).^(ii+1)'.*(ii'.^2 + ii'.^4).*U(ii+1))/N + .5*(-1)^(N+1)*N*U(N+1);
end