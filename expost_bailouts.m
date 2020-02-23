function [Ev_w_eb, Ev_b_eb] = expost_bailouts(frac, thresh, c_v_w, c_v_b, eq, param, glob, options)
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
ratiop          = glob.s(:, 2);

% [c_v_w, c_v_b]      = compute_value_function(eq, param, glob, options);
v_w                 = glob.Phi * c_v_w;
v_b                 = glob.Phi * c_v_b;

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder(1), glob.s(:, 1));           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder(2), ...
            ratiop .* (1 - frac .* (glob.s(:, 2) > thresh & glob.s(:, 3) < 1)));       % Basis for all b'
Phi_KBZp        = dprod(glob.Phi_Z, dprod(Phi_Bp, Phi_Kp)); 

v_w_eb              = Phi_KBZp * c_v_w;
v_b_eb              = Phi_KBZp * c_v_b;
Ev_w_eb             = utility_inv(sum(eq.dist .* v_w_eb), 0, param, glob, options);
Ev_b_eb             = utility_inv(sum(eq.dist .* v_b_eb), 0, param, glob, options);


end