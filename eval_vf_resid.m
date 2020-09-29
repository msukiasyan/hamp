function [res, jac, eq]  = eval_vf_resid(c, cnum, param, glob, options)
ns              = size(glob.s, 1);
cold            = c(1:ns);
ecold           = c(ns+1:2*ns);

%% Compute the optimal policies
optim_fast              = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'trust-region', ... %levenberg-marquardt
    'SpecifyObjectiveGradient',true, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-30, ...
    'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIterations', 2000);
cnumnew                 = fsolve(@(x) eval_vf_resid_backward_FR(ecold, x, param, glob, options), cnum, optim_fast);
[~, ~, eq]              = eval_vf_resid_backward_FR(ecold, cnumnew, param, glob, options);
%% Unpack
Phi             = glob.Phi;
Emat            = glob.Emat;
Kp              = eq.Kp;
Bp              = eq.Bp;
mup             = eq.mup;
c_w             = eq.c_w;
c_b             = eq.c_b;
L               = eq.L;
mu              = glob.s(:, 3);
B               = glob.s(:, 1) .* glob.s(:, 2);

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
% next_basis      = Phi_KBmuZp * PhiinvEmat;

%% Compute values from the functional equation
res(1:ns)       = Phi * cold - (glob.lambda * u_b + (1 - glob.lambda) * u_w - ((-mu) - (-mup)) .* B .* u_w_c + (-mup) .* (-c_w .* u_w_c - L .* u_w_l) ...
            + glob.beta_w * Phi_KBmuZp * ecold);
res(ns+1:2*ns)  = Phi * ecold - Emat * cold;
%% The Jacobian follows immediately by the envelope theorem
jac             = [Phi,     -glob.beta_w * Phi_KBmuZp;
                   -Emat,                       Phi];

end