function v_w = worker_welfare_pol(target_lev, slope, init_state, param, glob, options)
fprintf('Target lev: %3.4f\tSlope: %3.4f\n', target_lev, slope);

%% Compute
glob.target_lev     = target_lev;
glob.slope_pol      = slope;
glob.transf         = glob.slope_pol * (glob.s(:, 2) >= glob.target_lev) .* (glob.s(:, 3) < 2) .* abs(glob.s(:, 2) - glob.target_lev);
glob.bound_pol      = Inf;
options.save_eqbm   = 'N';
options.optim       = optimoptions(options.optim, 'Display', 'off');
eqbm_pol            = solve_eqbm(param, glob, options);

% Compute welfare
[c_v_w_pol, ~]      = compute_value_function(eqbm_pol, param, glob, options);

v_w                 = funeval(c_v_w_pol, glob.fspace, init_state);

end