function out = production_kl(z, k, l, param, glob, options)
    out         = param.alpha * (1 - param.alpha) * z .* ((l ./ k) .^ ( - param.alpha)) ./ k;
end