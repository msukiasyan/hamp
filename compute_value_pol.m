function [c_v_w_pl, c_v_b_pl] = compute_value_pol(eq, param, glob, options)
%% Unpack
Phi         = glob.Phi;
Emat        = glob.Emat;
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
jac             = [Phi,     -glob.beta_w * Phi_KBmuZp;
                   -Emat,                       Phi];
c_v_w_pl        = jac \ [u_w; zeros(size(glob.s, 1), 1)];
c_v_w_pl        = c_v_w_pl(1:glob.Ns);
c_v_b_pl        = jac \ [u_b; zeros(size(glob.s, 1), 1)];
c_v_b_pl        = c_v_b_pl(1:glob.Ns);

end