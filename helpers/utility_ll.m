function out = utility_ll(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = - param.chi * param.nu * l .^ (param.nu - 1);
    else
        out         = - param.chi * param.nu * l .^ (param.nu - 1) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma) ...
            + (param.chi * l .^ (param.nu)) .^ 2 .* (-param.gamma) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 1);
    end
end