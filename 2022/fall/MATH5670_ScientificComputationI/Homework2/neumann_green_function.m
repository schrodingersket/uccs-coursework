function G = neumann_green_function(xbar, x)
    arguments
        xbar   (1,1) double {mustBeNonnegative,mustBeLessThanOrEqual(xbar, 1)}
        x      (1,:) double {mustBeNonnegative,mustBeLessThanOrEqual(x, 1)}
    end

    G = zeros(1, length(x));
    for i=1:1:length(x)
        if x(i) <= xbar
            G(i) = (xbar - 1);
        else
            G(i) = x(i) - 1;
    end
end
