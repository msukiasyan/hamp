function out = cap_prod_prime_prime(x, param, glob, options)
    out         = glob.A * (glob.xi - 1) * (glob.A * (x - param.delta) + 1) .^ (glob.xi - 2);
end