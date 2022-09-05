function G = green_function(xbar, x, f_xbar)
    arguments
        xbar   (1,1) double {mustBeNonnegative,mustBeLessThanOrEqual(x_bar, 1)}
        x      (1,:) double {mustBeNonnegative,mustBeLessThanOrEqual(x, 1)}
        f_xbar (1,1) double
    end

    gf = zeros(1, length(x));
    for i=1:1:length(x)
        if x(i) <= xbar
            gf(i) = (xbar - 1) .* x(i);
        else
            gf(i) = xbar .* (x(i) - 1);
    end

    G = f_xbar .* gf;
end
