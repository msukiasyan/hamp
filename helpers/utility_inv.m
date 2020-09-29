function c = utility_inv(u, l, param, glob, options)
    if options.GHH == 'N'
        c           = ((1 - param.gamma) * u) .^ (1 / (1 - param.gamma));
    else
        c           = ((1 - param.gamma) * u) .^ (1 / (1 - param.gamma)) + param.chi * l .^ (1 + param.nu) / (1 + param.nu);
    end
end