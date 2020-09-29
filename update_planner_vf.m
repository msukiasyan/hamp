function [c_v] = update_planner_vf(eq, EWc, param, glob, options)
%% Unpack
Phi         = glob.Phi;
% Phiinv      = glob.Phiinv;
Emat        = glob.Emat;
% PhiinvEmat  = glob.PhiinvEmat;
Kp          = eq.Kp;
Bp          = eq.Bp;
mup         = eq.mup;
c_w         = eq.c_w;
c_b         = eq.c_b;
L           = eq.L;
mu          = glob.s(:, 3);
B           = glob.s(:, 1) .* glob.s(:, 2);

%% Prepare ingredients
ratiop          = Bp ./ Kp;

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop);       % Basis for all b'
Phi_mup         = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup);       % Basis for all mu'
Phi_KBmuZp      = dprod(glob.Phi_Z, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp)));                     % Basis for all (k', b', mu', z')
u_w             = utility(c_w, L, param, glob, options);
u_w_c           = utility_c(c_w, L, param, glob, options);
u_w_l           = utility_l(c_w, L, param, glob, options);
u_b             = utility(c_b, 0, param, glob, options);

%% Compute values from the functional equation
c_v             = Phi \ (glob.beta_w * Phi_KBmuZp * EWc + ...
    (glob.lambda * u_b + (1 - glob.lambda) * u_w - ((-mu) - (-mup)) .* B .* u_w_c + (-mup) .* (-c_w .* u_w_c - L .* u_w_l)));

end