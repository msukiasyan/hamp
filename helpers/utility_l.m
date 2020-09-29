function out = utility_l(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = - param.chi * l .^ (param.nu);
    else
        out         = - param.chi * l .^ (param.nu) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma);
    end
end