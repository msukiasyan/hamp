function out = production_l(z, k, l, param, glob, options)
    out         = (1 - param.alpha) * z .* (k ./ l) .^ (param.alpha);
end