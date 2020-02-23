function out = cap_prod(x, param, glob, options)
    out         = ((glob.A * (x - param.delta) + 1) .^ glob.xi - 1) / glob.A / glob.xi + param.delta;
end