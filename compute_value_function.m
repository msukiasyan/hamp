function [c_v_w, c_v_b] = compute_value_function(eq, param, glob, options, perc_c_w, perc_c_b)
%% Unpack
Phi         = glob.Phi;
Phiinv      = glob.Phiinv;
Emat        = glob.Emat;
Kp          = eq.Kp;
Bp          = eq.Bp;
c_w         = eq.c_w;
c_b         = eq.c_b;
L           = eq.L;

%% Prepare ingredients
ratiop          = Bp ./ Kp;

if nargin <= 4
    perc_c_w    = 1;
    perc_c_b    = 1;
end

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder(2), ratiop);       % Basis for all b'
Phi_KBZp        = dprod(glob.Phi_Z, dprod(Phi_Bp, Phi_Kp));                 % Basis for all (k', b', z')
u_w             = utility(c_w * perc_c_w, L, param, glob, options);
u_b             = utility(c_b * perc_c_b, 0, param, glob, options);

%% Compute values from the functional equation
c_v_w           = (Phi - glob.beta_w * Phi_KBZp * Phiinv * Emat) \ u_w;
c_v_b           = (Phi - glob.beta_b * Phi_KBZp * Phiinv * Emat) \ u_b;

end