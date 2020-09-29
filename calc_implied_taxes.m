function eq     = calc_implied_taxes(eq, param, glob, options)
%% Unpack
Phi         = glob.Phi;
Emat        = glob.Emat;
Kp          = eq.Kp;
Bp          = eq.Bp;
mup         = eq.mup;
c_w         = eq.c_w;
c_b         = eq.c_b;
L           = eq.L;
K           = glob.s(:, 1);
B           = glob.s(:, 1) .* glob.s(:, 2);
mu          = glob.s(:, 3);
Z           = glob.s(:, 4);
ns          = size(glob.s, 1);
%% Prepare objects
ratiop          = Bp ./ Kp;

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop);       % Basis for all b'
Phi_mup         = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup);       % Basis for all mu'
Phi_KBmuZp      = dprod(glob.Phi_Z, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp)));                     % Basis for all (k', b', mu', z')

u_w_c           = utility_c(c_w, L, param, glob, options);
u_b_c           = utility_c(c_b, zeros(ns, 1), param, glob, options);
u_w_l           = utility_l(c_w, L, param, glob, options);
%% Compute implied tax rates
eq.implied_tax_lab  = 1 + u_w_l ./ u_w_c ./ production_l(Z, K, L, param, glob, options);

PhiEmu_w            = Phi \ (Emat * (Phi \ u_w_c));
PhiEmu_b            = Phi \ (Emat * (Phi \ u_b_c));
Emu_wp              = Phi_KBmuZp * PhiEmu_w;
Emu_bp              = Phi_KBmuZp * PhiEmu_b;
implied_r           = 1 ./ (glob.beta .* Emu_wp ./ u_w_c);
implied_r_b         = 1 ./ (glob.beta .* Emu_bp ./ u_b_c);

eq.implied_tax_r    = 1 - implied_r ./ implied_r_b;
eq.transfer         = eq.implied_tax_lab .* L + eq.implied_tax_r .* B ./ implied_r_b;
eq.r                = implied_r;

end