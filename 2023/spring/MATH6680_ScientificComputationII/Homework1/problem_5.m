set(0,'DefaultFigureWindowStyle','docked')
% p3.m - band-limited interpolation
hh = [2^-3, 2^-4, 2^-5, 2^-6];
max_errs = zeros(2, length(hh));

square_wave = @(t) (abs(t) <= 3);
hat_function = @(t) max(0, 1 - abs(t)/3);
clf('reset');
for k = 1:length(hh)
%     figure;
    h = hh(k);
    xmax = 10;
    x = -xmax:h:xmax; % computational grid
    xx = -xmax-h/20:h/10:xmax+h/20; % plotting grid
    for plt = 1:2
%         subplot(2, 1, plt)
        switch plt
            case 1
                v = square_wave(x);
                vv = square_wave(xx);
            case 2
                v = hat_function(x);
                vv = hat_function(xx);
        end

        p = zeros(size(xx));
        for i = 1:length(x)
            p = p + v(i)*sin(pi*(xx-x(i))/h)./(pi*(xx-x(i))/h);
        end
%         
%         hold on;
%         plot(xx, p, 'b');
%         plot(x, v, 'k');

        max_errs(plt, k) = norm(p - vv, Inf);
    end
end

hh
max_errs

subplot(2, 1, 1);
loglog(hh, max_errs(1, :));
title('Square Wave Error')
xlabel('h')
ylabel('error')


subplot(2, 1, 2);
loglog(hh, max_errs(2, :));
title('Hat Function Error')
xlabel('h')
ylabel('error')

print('-dpng', 'problem_5')
