function out = utility_cc(c, l, param, glob, options)
    out         = (- param.gamma) * c .^ (- param.gamma - 1);
end