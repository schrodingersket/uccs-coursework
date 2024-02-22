% p5.m - repetition of p4.m via FFT
% For complex v, delete "real" commands.
% Differentiation of a hat function:

v = @(x) exp(sin(x));

% First, time FFT for powers of 2
%
K = (2 * ones(1, 16)).^(0:15);
TK = zeros(size(K));

for i=1:length(K)
    N = K(i);
    h = 2*pi/N;
    TK(i) = time_fft(v, h, N);
end

subplot(2, 1, 1);
plot(K, TK, '-k');
title('FFT for N=2^k');
xlabel('N');
ylabel('Time (s)');

NN = 500:520;
TN = zeros(size(NN));

for i=1:length(NN)
    N = NN(i);
    h = 2*pi/N;
    TN(i) = time_fft(v, h, N);
end

subplot(2, 1, 2);
plot(NN, TN, '-k');
title('FFT for N~500');
xlabel('N');
ylabel('Time (s)');
print('-dpng', 'problem_6');


function elapsed = time_fft(v, h, N)
    xx = h*(1:N)';
    vv = v(xx);
    iterations = 100;
    tic
    for i=1:iterations
        fft(vv);
    end
    elapsed = toc / iterations;
end