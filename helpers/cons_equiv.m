function [ceq_w, ceq_b] = cons_equiv(v_w, v_b, eq, param, glob, options)
%% Unpack
Phi         = glob.Phi;
Phiinv      = glob.Phiinv;
Emat        = glob.Emat;
Kp          = eq.Kp;
Bp          = eq.Bp;
c_w         = eq.c_w;
c_b         = eq.c_b;
L           = eq.L;

%% Compute
Np                  = 50;
perc_grid           = linspace(0.95, 1.15, Np);
vals_w              = zeros(Np, glob.Ns);
vals_b              = zeros(Np, glob.Ns);

for pi = 1:Np
    [c_v_w, c_v_b]  = compute_value_function(eq, param, glob, options, perc_grid(pi), perc_grid(pi));
    vals_w(pi, :)   = glob.Phi * c_v_w;
    vals_b(pi, :)   = glob.Phi * c_v_b;
end

f                   = @(x) arrayfun(@(y) interp1(perc_grid, vals_w(:, y), x(y)), 1:glob.Ns)' - v_w;
ceq_w               = binsearchx(f, ones(glob.Ns, 1) * 0.95, ones(glob.Ns, 1) * 1.15);

f                   = @(x) arrayfun(@(y) interp1(perc_grid, vals_b(:, y), x(y)), 1:glob.Ns)' - v_b;
ceq_b               = binsearchx(f, ones(glob.Ns, 1) * 0.95, ones(glob.Ns, 1) * 1.15);

end