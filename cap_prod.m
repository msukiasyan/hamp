function out = cap_prod(x, param, glob, options)
    out         = (glob.A * (x - param.delta) + 1) .^ glob.xi - 1 + param.delta;
end