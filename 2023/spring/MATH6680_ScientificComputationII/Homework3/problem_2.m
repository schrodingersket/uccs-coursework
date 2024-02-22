p = 2;
z = @(theta) exp(1i*theta);
f = @(z) log(1 + z./p);

M = 20;  % Number of Taylor series coefficients
maxN = 2*M; % Number of points to use for FFT

% Analytic computation of Taylor coefficients of f
%
taylor_coefficients = [0 -1*(-1/p).^(1:M)./(1:M) ];

err = zeros(1, maxN);

for k=2:maxN
    tt = 2*pi*(0:k-1)/k;
    zz = z(tt);
    
    taylor_fft = real(fft(f(zz), k) / length(zz));
    err(k) = abs(taylor_coefficients(2) - taylor_fft(2));
end

semilogy(2:maxN, err(2:end), '-o')
title('Error in Taylor coefficient a_1')
xlabel('N')
ylabel('Error')

print('problem_2.png', '-dpng')

disp('Ratio of sequential terms (convergence factor)')
fprintf('%f \n', err(2:end)./err(1:end-1))
disp('')

% FFT computation
%
tt = 2*pi*(0:maxN-1)/maxN;
zz = z(tt);
taylor_fft = real(fft(f(zz), maxN) / length(zz));

fprintf('Computation via FFT with N=%d: \n', maxN)
for k = 1:M
    fprintf('a_%d: \t %d \n', k-1, taylor_fft(k));
end