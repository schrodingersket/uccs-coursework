function out=bestfitnormal(data)
    [M, dim]=size(data);

    Const = ones(M,1); % Vector of the constant term in the RHS
    out = data\Const; % Find the coefficients
end
