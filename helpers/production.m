function out = production(z, k, l, param, glob, options)
    out         = z .* k .^ param.alpha .* l .^ (1 - param.alpha);
end