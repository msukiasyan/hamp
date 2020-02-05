function out = cap_prod_prime(x, param, glob, options)
    out         = glob.A * glob.xi * (glob.A * (x - param.delta) + 1) .^ (glob.xi - 1);
end