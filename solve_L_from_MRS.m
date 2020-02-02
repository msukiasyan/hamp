function out = solve_L_from_MRS(z, k, c_w, param, glob, options)
    out         = ((1 - param.alpha) / param.chi * z .* k .* param.alpha .* c_w .^ (-param.gamma)) .^ (1 / (param.alpha + param.nu));
end