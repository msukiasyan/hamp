function out = utility_ccc(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = (- param.gamma) * (- param.gamma - 1) * c .^ (- param.gamma - 2);
    else
        out         = (- param.gamma) * (- param.gamma - 1) * (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 2);
    end
end