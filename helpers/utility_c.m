function out = utility_c(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = c .^ (- param.gamma);
    else
        out         = (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma);
    end
end