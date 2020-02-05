function out = production_k(z, k, l, param, glob, options)
    out         = param.alpha * z .* (l ./ k) .^ (1 - param.alpha);
end