nu = 4;
h = 0.1;

% starting and ending time
tstart = 0;
tstop = 2 / h * pi;


% make the time points
t = linspace(tstart, tstop);

% calculate z
numerator = @(xi) 1 - nu/2 * (1 - exp(-i * xi * h));
denominator = @(xi) 1 + nu/2 * (1 - exp(-i * xi * h));

z = numerator(t) ./ denominator(t);

% extract real and imaginary parts
rr = real(z);
ii = imag(z);

scatter(rr, ii);
axis([-1 1 -1 1]);

print('-dpng', 'problem_2ci')
