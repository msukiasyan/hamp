function out = utility_ccl(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = zeros(size(c, 1), 1);
    else
        out         = (- param.gamma) * (- param.gamma - 1) * (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 2) .* ...
            (-param.chi * l .^ param.nu);
    end
end