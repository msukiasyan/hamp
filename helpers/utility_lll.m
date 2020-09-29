function out = utility_lll(c, l, param, glob, options)
    if options.GHH == 'N'
        out         = - param.chi * (param.nu - 1) * param.nu * l .^ (param.nu - 2);
    else
        out         = - param.chi * param.nu * (param.nu - 1) * l .^ (param.nu - 2) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma) + ...
            param.gamma * (- param.chi) * param.nu * l .^ (param.nu - 1) * (- param.chi) * l .^ (param.nu) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 1) + ...
            2 * param.chi * param.nu * l .^ (param.nu - 1) .* (param.chi * l .^ (param.nu)) .* (-param.gamma) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 1) + ...
            (-param.chi * l .^ (param.nu)) .^ 3 .* (-param.gamma) * (-param.gamma - 1) .* (c - param.chi * l .^ (1 + param.nu) / (1 + param.nu)) .^ (- param.gamma - 2);
    end
end