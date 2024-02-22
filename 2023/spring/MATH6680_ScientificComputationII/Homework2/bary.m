function b = bary(x, u, xx)
    top = zeros(size(xx));
    bottom = zeros(size(xx));
    for k = 1:(length(x))
        a_k = prod(x(k) - x([1:(k-1), (k+1):end]));
        common = 1./((xx - x(k)).*a_k);
        top = top + common * u(k);
        bottom = bottom + common;
    end

    b = top ./ bottom;
end
