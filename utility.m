function out = utility(c, l, param, glob, options)
    out         = c .^ (1 - param.gamma) / (1 - param.gamma) - param.chi * l .^ (1 + param.nu) / (1 + param.nu);
end