function [vs] = planner_vf(c_0, wc, param, glob, options, s)
%% Unpack
Phi         = glob.Phi;
Phiinv      = glob.Phiinv;
Emat        = glob.Emat;

if nargin < 6
    s       = glob.s;
end

ns          = size(s, 1);

c_w         = c_0(1:ns);
L           = c_0(ns + 1:2 * ns);
Bp          = c_0(2 * ns + 1:3 * ns);
mu          = s(:, 3);
B           = s(:, 1) .* s(:, 2);
Z           = s(:, 4);
K           = s(:, 1);

H_l             = (B - c_w) .* utility_cl(c_w, L, param, glob, options) -  ...
    utility_ll(c_w, L, param, glob, options) .* L - utility_l(c_w, L, param, glob, options);
H_c             = (B - c_w) .* utility_cc(c_w, L, param, glob, options) -  ...
    utility_c(c_w, L, param, glob, options) - utility_cl(c_w, L, param, glob, options) .* L;

gam             = -((production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_l(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cl(c_w, L, param, glob, options)) ./ (-H_l - production_l(Z, K, L, param, glob, options) .* H_c));
c_b             = utility_c_inv(max(((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + gam .* H_c) / glob.lambda, 0), 0, param, glob, options);
Y               = production(Z, K, L, param, glob, options);
I               = Y - c_b - c_w;
q               = 1 ./ cap_prod_prime(I ./ K, param, glob, options);
Kp              = K .* (cap_prod(I ./ K, param, glob, options) + 1 - param.delta);
mup             = gam;
ratiop          = Bp ./ Kp;


%% Prepare ingredients
Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop);       % Basis for all b'
Phi_mup         = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup);       % Basis for all mu'
Phi_Z           = splibas(glob.zgrid0, 0, glob.spliorder_FR(4), Z);       % Basis for all mu'
Phi_KBmuZp      = dprod(Phi_Z, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp)));                     % Basis for all (k', b', mu', z')
u_w             = utility(c_w, L, param, glob, options);
u_w_c           = utility_c(c_w, L, param, glob, options);
u_w_l           = utility_l(c_w, L, param, glob, options);
u_b             = utility(c_b, 0, param, glob, options);

%% Compute values from the functional equation
vs             = glob.beta_w * Phi_KBmuZp * Phiinv * Emat * wc + ...
    (glob.lambda * u_b + (1 - glob.lambda) * u_w - (mu - mup) .* B .* u_w_c + mup .* (-c_w .* u_w_c - L .* u_w_l));

end