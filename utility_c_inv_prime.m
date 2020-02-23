function out = utility_c_inv_prime(mu, l, param, glob, options)
    if options.GHH == 'N'
        out         = (-1 / param.gamma) * mu .^ (-1 / param.gamma - 1);
    else
        out         = (-1 / param.gamma) * mu .^ (-1 / param.gamma - 1);
    end
end