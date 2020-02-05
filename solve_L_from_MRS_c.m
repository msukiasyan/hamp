function out = solve_L_from_MRS_c(z, k, c_w, param, glob, options)
    out         = ((1 - param.alpha) / param.chi * z .* k .* param.alpha) .^ (1 / (param.alpha + param.nu)) .* ...
                    (-param.gamma / (param.alpha + param.nu)) .* c_w .^ (-param.gamma / (param.alpha + param.nu) - 1);
end