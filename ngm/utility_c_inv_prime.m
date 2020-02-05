function out = utility_c_inv_prime(mu, l, param, glob, options)
    out         = (-1 / param.gamma) * mu .^ (-1 / param.gamma - 1);
end