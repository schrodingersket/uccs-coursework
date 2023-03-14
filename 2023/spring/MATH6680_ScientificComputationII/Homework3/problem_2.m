z = @(theta) exp(1i*theta);
f = @(z) log(1 + z./2);

N = 50;
M = 20;  % Number of Taylor series coefficients

% Manual computation via periodic trapezoid rule
%
taylor_trapsum = zeros(1, M);

tt = 2*pi*(0:N-1)/N;
zz = z(tt);

for m = 0:M-1
    taylor_trapsum(m+1) = real(mean(exp(-1i * m * tt) .* f(zz)));
end

disp('Computation via periodic trapezoidal sum:')
disp(taylor_trapsum(1:5))

% FFT computation
%
taylor_fft = real(fft(f(zz), N) / length(zz));

disp('Computation via FFT:')
disp(taylor_fft(1:5))

% Print sequential ratios
%
% taylor_fft(2:end) ./ taylor_fft(1:end-1)
