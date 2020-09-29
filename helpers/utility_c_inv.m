function out = utility_c_inv(mu, l, param, glob, options)
    if options.GHH == 'N'
        out         = mu .^ (-1 / param.gamma);
    else
        out         = mu .^ (-1 / param.gamma) + param.chi * l .^ (1 + param.nu) / (1 + param.nu);
    end
end