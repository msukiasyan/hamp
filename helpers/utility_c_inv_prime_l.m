function out = utility_c_inv_prime_l(mu, l, param, glob, options)
    if options.GHH == 'N'
        out         = 0;
    else
        out         = param.chi * l .^ param.nu;
    end
end