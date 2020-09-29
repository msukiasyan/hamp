function out = utility_cl(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = zeros(size(c, 1), 1);
    else
        out         = (- param.gamma) * (-param.chi * l .^ param.nu) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 1);
    end
end