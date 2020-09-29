function out = production_ll(z, k, l, param, glob, options)
    out         = (1 - param.alpha) * (-param.alpha) * z .* k .^ (param.alpha) .* l .^ (-param.alpha - 1);
end