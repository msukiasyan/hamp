function out = utility_cc(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = (- param.gamma) * c .^ (- param.gamma - 1);
    else
        out         = (- param.gamma) * (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 1);
    end
end