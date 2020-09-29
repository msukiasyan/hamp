function out = solve_L_from_MRS(z, k, c_w, param, glob, options)
    if options.GHH == 'N'
        out         = ((1 - param.alpha) / param.chi * z .* k .^ param.alpha .* c_w .^ (-param.gamma)) .^ (1 / (param.alpha + param.nu));
    else
        out         = ((1 - param.alpha) / param.chi * z .* k .^ param.alpha) .^ (1 / (param.alpha + param.nu));
    end
end