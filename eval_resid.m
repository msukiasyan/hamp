function [res, jac, L, Y, I, q, Kp, Bp, r] = eval_resid(c, param, glob, options)
%% Globals 
s               = glob.s;  
ns              = size(s, 1);
Phi             = glob.Phisp;
Phiinv          = glob.Phiinv;
Phiu            = glob.Phiu;
Phil            = glob.Phil;
basiscast       = glob.basiscast;
K               = glob.s(:, 1);
B               = glob.s(:, 1) .* glob.s(:, 2);
Z               = glob.s(:, 3);
% Unpack
c1              = c(1:ns);
c2              = c(ns + 1:2 * ns);
c3              = c(2 * ns + 1:end);
c_w             = Phi * c1;
c_b             = Phi * c2;
r               = Phi * c3;

% c_w(c_w < 1e-5)    = 1e-5;
% c_b(c_b < 1e-5)    = 1e-5;

%% Solve equations 
L               = solve_L_from_MRS(Z, K, c_w, param, glob, options);
Y               = production(Z, K, L, param, glob, options);
I               = Y - c_b - c_w;
q               = 1 ./ cap_prod_prime(I ./ K, param, glob, options);
Kp              = K .* (cap_prod(I ./ K, param, glob, options) + 1 - param.delta);
mpk             = production_k(Z, K, L, param, glob, options);
Bp              = r .* (B + c_b + q .* Kp - mpk .* K - (1 - param.delta) * q .* K);
mu_w            = utility_c(c_w, L, param, glob, options);
mu_b            = utility_c(c_b, L, param, glob, options);
mu_bprod        = mu_b .* (mpk + (1 - param.delta) * q);

% Enforce bounds
% Kp              = min(max(Kp, glob. kmin), glob.kmax);
% Bp              = min(max(Bp, glob. bmin), glob.bmax);

%% Create basis matrices for next states
ratiop          = Bp ./ Kp;

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder(2), ratiop);           % Basis for all b'
Phi_KBZp        = dprod(glob.Phi_Z, dprod(Phi_Bp, Phi_Kp));                 % Basis for all (k', b', z')

%% Compute next state objects
PhiEmu_bprod    = basiscast * mu_bprod;                                     % Approximate on the same basis, take expectation and get coefficients again
PhiEmu_b        = basiscast * mu_b;                     
PhiEmu_w        = basiscast * mu_w;                     
Emu_bprodp      = Phi_KBZp * PhiEmu_bprod;                               % Evaluate at (k', b', z)
Emu_bp          = Phi_KBZp * PhiEmu_b;                                   % Evaluate at (k', b', z)
Emu_wp          = Phi_KBZp * PhiEmu_w;                                   % Evaluate at (k', b', z)


%% Compute residuals
res                 = zeros(3 * ns, 1);
res(1:ns)           = (c_w) - utility_c_inv(glob.beta_w * r .* Emu_wp, L, param, glob, options);                     % Euler equation for workers
res(ns+1:2*ns)      = (c_b) - utility_c_inv(glob.beta_b * r .* Emu_bp, zeros(ns, 1), param, glob, options);                     % Euler equation for bankers
res(2*ns+1:3*ns)    = utility_c_inv(Emu_bprodp, zeros(ns, 1), param, glob, options) - ...
                            utility_c_inv(Emu_bp .* q .* r, zeros(ns, 1), param, glob, options);                        % Arbitrage

%% Compute the jacobian if requested
jac             = [];
if (nargout == 2)
    % Compute derivatives with respect to all three inputs
    L_c_w           = solve_L_from_MRS_c(Z, K, c_w, param, glob, options);
    L_c_b           = zeros(ns, 1);
    L_r             = zeros(ns, 1);
    
    Y_c_w           = production_l(Z, K, L, param, glob, options) .* L_c_w; % Chain rule
    Y_c_b           = zeros(ns, 1);
    Y_r             = zeros(ns, 1);
    
    I_c_w           = Y_c_w - ones(ns, 1);
    I_c_b           = Y_c_b - ones(ns, 1);
    I_r             = zeros(ns, 1);
    
    q_c_w           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c_w;
    q_c_b           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c_b;
    q_r             = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_r;
    
    Kp_c_w          = (1 ./ q) .* I_c_w;
    Kp_c_b          = (1 ./ q) .* I_c_b;
    Kp_r            = (1 ./ q) .* I_r;
    
    mpk_c_w         = production_kl(Z, K, L, param, glob, options) .* L_c_w;
    mpk_c_b         = production_kl(Z, K, L, param, glob, options) .* L_c_b;
    mpk_r           = production_kl(Z, K, L, param, glob, options) .* L_r;
    
    Bp_c_w          = r .* (q_c_w .* Kp + q .* Kp_c_w - mpk_c_w .* K - (1 - param.delta) * q_c_w .* K);
    Bp_c_b          = r .* (ones(ns, 1) + q_c_b .* Kp + q .* Kp_c_b - mpk_c_b .* K - (1 - param.delta) * q_c_b .* K);
    Bp_r            = ones(ns, 1) .* (B + c_b + q .* Kp - mpk .* K - (1 - param.delta) * q .* K) + ...
                        r .* (q_r .* Kp + q .* Kp_r - mpk_r .* K - (1 - param.delta) * q_r .* K);
                    
    %-------------
    Bp_c_w          = Bp_c_w ./ Kp - Bp ./ (Kp .^2) .* Kp_c_w;
    Bp_c_b          = Bp_c_b ./ Kp - Bp ./ (Kp .^2) .* Kp_c_b;
    Bp_r            = Bp_r ./ Kp - Bp ./ (Kp .^2) .* Kp_r;
    %-------------
    
    mu_w_c_w        = utility_cc(c_w, L, param, glob, options) .* ones(ns, 1) + utility_cl(c_w, L, param, glob, options) .* L_c_w;
    mu_w_c_b        = utility_cl(c_w, L, param, glob, options) .* L_c_b;
    mu_w_r          = utility_cl(c_w, L, param, glob, options) .* L_r;
    
    mu_b_c_w        = utility_cl(c_b, L, param, glob, options) .* L_c_w;
    mu_b_c_b        = utility_cc(c_b, L, param, glob, options) .* ones(ns, 1) + utility_cl(c_b, L, param, glob, options) .* L_c_b;
    mu_b_r          = utility_cl(c_b, L, param, glob, options) .* L_r;
    
    mu_bprod_c_w    = mu_b_c_w .* (mpk + (1 - param.delta) * q) + mu_b .* (mpk_c_w + (1 - param.delta) * q_c_w);
    mu_bprod_c_b    = mu_b_c_b .* (mpk + (1 - param.delta) * q) + mu_b .* (mpk_c_b + (1 - param.delta) * q_c_b);
    mu_bprod_r      = mu_b_r .* (mpk + (1 - param.delta) * q) + mu_b .* (mpk_r + (1 - param.delta) * q_r);
    
    % Create derivative matrices
    Phi_Kp_der      = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp, 1);
    Phi_Bp_der      = splibas(glob.bgrid0, 0, glob.spliorder(2), ratiop, 1);
    Phi_KBZp_Kder   = dprod(glob.Phi_Z, dprod(Phi_Bp, Phi_Kp_der));      
    Phi_KBZp_Bder   = dprod(glob.Phi_Z, dprod(Phi_Bp_der, Phi_Kp));  
    
    Kder_bprod      = spdiags(Phi_KBZp_Kder * PhiEmu_bprod, 0, ns, ns);
    Bder_bprod      = spdiags(Phi_KBZp_Bder * PhiEmu_bprod, 0, ns, ns);
    Kder_b          = spdiags(Phi_KBZp_Kder * PhiEmu_b, 0, ns, ns);
    Bder_b          = spdiags(Phi_KBZp_Bder * PhiEmu_b, 0, ns, ns);
    Kder_w          = spdiags(Phi_KBZp_Kder * PhiEmu_w, 0, ns, ns);
    Bder_w          = spdiags(Phi_KBZp_Bder * PhiEmu_w, 0, ns, ns);
    
    diagKp_c_w      = spdiags(Kp_c_w, 0, ns, ns);
    diagKp_c_b      = spdiags(Kp_c_b, 0, ns, ns);
    diagKp_r        = spdiags(Kp_r, 0, ns, ns);
    diagBp_c_w      = spdiags(Bp_c_w, 0, ns, ns);
    diagBp_c_b      = spdiags(Bp_c_b, 0, ns, ns);
    diagBp_r        = spdiags(Bp_r, 0, ns, ns);
    
    next_basis      = Phi_KBZp * basiscast;
   

    Emu_bprodp_c_w  = Kder_bprod * diagKp_c_w + Bder_bprod * diagBp_c_w + ... 
                        next_basis * spdiags(mu_bprod_c_w, 0, ns, ns);
    Emu_bprodp_c_b  = Kder_bprod * diagKp_c_b + Bder_bprod * diagBp_c_b + ...
                        next_basis * spdiags(mu_bprod_c_b, 0, ns, ns);
    Emu_bprodp_r    = Kder_bprod * diagKp_r + Bder_bprod * diagBp_r + ...
                        next_basis * spdiags(mu_bprod_r, 0, ns, ns);
                    
    Emu_bp_c_w      = Kder_b * diagKp_c_w + Bder_b * diagBp_c_w + ...
                        next_basis * spdiags(mu_b_c_w, 0, ns, ns);
    Emu_bp_c_b      = Kder_b * diagKp_c_b + Bder_b * diagBp_c_b + ...
                        next_basis * spdiags(mu_b_c_b, 0, ns, ns);
    Emu_bp_r        = Kder_b * diagKp_r + Bder_b * diagBp_r + ...
                        next_basis * spdiags(mu_b_r, 0, ns, ns);
    
    Emu_wp_c_w      = Kder_w * diagKp_c_w + Bder_w * diagBp_c_w + ...
                        next_basis * spdiags(mu_w_c_w, 0, ns, ns);
    Emu_wp_c_b      = Kder_w * diagKp_c_b + Bder_w * diagBp_c_b + ...
                        next_basis * spdiags(mu_w_c_b, 0, ns, ns);
    Emu_wp_r        = Kder_w * diagKp_r + Bder_w * diagBp_r + ...
                        next_basis * spdiags(mu_w_r, 0, ns, ns);

    % Calculate the blocks of the jacobian
%     res1_c_w        = spdiags(mu_w_c_w, 0, ns, ns) - glob.beta_w * ...
%                     (spdiags(r, 0, ns, ns) * Emu_wp_c_w);
%     res1_c_b        = spdiags(mu_w_c_b, 0, ns, ns) - glob.beta_w * ...
%                     (spdiags(r, 0, ns, ns) * Emu_wp_c_b);
%     res1_r          = spdiags(mu_w_r, 0, ns, ns) - glob.beta_w * ...
%                     (spdiags(r, 0, ns, ns) * Emu_wp_r + spdiags(Emu_wp, 0, ns, ns));
%                 
%     res2_c_w        = spdiags(mu_b_c_w, 0, ns, ns) - glob.beta_b * ...
%                     (spdiags(r, 0, ns, ns) * Emu_bp_c_w);
%     res2_c_b        = spdiags(mu_b_c_b, 0, ns, ns) - glob.beta_b * ...
%                     (spdiags(r, 0, ns, ns) * Emu_bp_c_b);
%     res2_r          = spdiags(mu_b_r, 0, ns, ns) - glob.beta_b * ...
%                     (spdiags(r, 0, ns, ns) * Emu_bp_r + spdiags(Emu_bp, 0, ns, ns));
%                 
%     res3_c_w        = Emu_bprodp_c_w - (spdiags(q .* r, 0, ns, ns) * Emu_bp_c_w + ...
%                         spdiags(Emu_bp .* r, 0, ns, ns) * spdiags(q_c_w, 0, ns, ns));
%     res3_c_b        = Emu_bprodp_c_b - (spdiags(q .* r, 0, ns, ns) * Emu_bp_c_b + ...
%                         spdiags(Emu_bp .* r, 0, ns, ns) * spdiags(q_c_b, 0, ns, ns));
%     res3_r          = Emu_bprodp_r - (spdiags(q .* r, 0, ns, ns) * Emu_bp_r + ...
%                         spdiags(Emu_bp .* r, 0, ns, ns) * spdiags(q_r, 0, ns, ns) + ...
%                         spdiags(Emu_bp .* q, 0, ns, ns));

    mu_prime1       = spdiags(utility_c_inv_prime(glob.beta_w * r .* Emu_wp, L, param, glob, options), 0, ns, ns);
    mu_prime2       = spdiags(utility_c_inv_prime(glob.beta_b * r .* Emu_bp, zeros(ns, 1), param, glob, options), 0, ns, ns);
    mu_prime31      = spdiags(utility_c_inv_prime(Emu_bprodp, zeros(ns, 1), param, glob, options), 0, ns, ns);
    mu_prime32      = spdiags(utility_c_inv_prime(Emu_bp .* q .* r, zeros(ns, 1), param, glob, options), 0, ns, ns);

    res1_c_w        = speye(ns, ns) - mu_prime1 * glob.beta_w * ...
                    (spdiags(r, 0, ns, ns) * Emu_wp_c_w);
    res1_c_b        =  - mu_prime1 * glob.beta_w * ...
                    (spdiags(r, 0, ns, ns) * Emu_wp_c_b);
    res1_r          =  - mu_prime1 * glob.beta_w * ...
                    (spdiags(r, 0, ns, ns) * Emu_wp_r + spdiags(Emu_wp, 0, ns, ns));
                
    res2_c_w        =  - mu_prime2 * glob.beta_b * ...
                    (spdiags(r, 0, ns, ns) * Emu_bp_c_w);
    res2_c_b        = speye(ns, ns) - mu_prime2 * glob.beta_b * ...
                    (spdiags(r, 0, ns, ns) * Emu_bp_c_b);
    res2_r          =  - mu_prime2 * glob.beta_b * ...
                    (spdiags(r, 0, ns, ns) * Emu_bp_r + spdiags(Emu_bp, 0, ns, ns));
                
    res3_c_w        = mu_prime31 * Emu_bprodp_c_w - mu_prime32 * (spdiags(q .* r, 0, ns, ns) * Emu_bp_c_w + ...
                        spdiags(Emu_bp .* r, 0, ns, ns) * spdiags(q_c_w, 0, ns, ns));
    res3_c_b        = mu_prime31 * Emu_bprodp_c_b - mu_prime32 * (spdiags(q .* r, 0, ns, ns) * Emu_bp_c_b + ...
                        spdiags(Emu_bp .* r, 0, ns, ns) * spdiags(q_c_b, 0, ns, ns));
    res3_r          = mu_prime31 * Emu_bprodp_r - mu_prime32 * (spdiags(q .* r, 0, ns, ns) * Emu_bp_r + ...
                        spdiags(Emu_bp .* r, 0, ns, ns) * spdiags(q_r, 0, ns, ns) + ...
                        spdiags(Emu_bp .* q, 0, ns, ns));
    
    % Put together the jacobian
    jac                             = zeros(3 * ns, 3 * ns);
    jac(1:ns, 1:ns)                 = res1_c_w * Phi;
    jac(1:ns, ns+1:2*ns)            = res1_c_b * Phi;
    jac(1:ns, 2*ns+1:3*ns)          = res1_r * Phi;
    jac(ns+1:2*ns, 1:ns)            = res2_c_w * Phi;
    jac(ns+1:2*ns, ns+1:2*ns)       = res2_c_b * Phi;
    jac(ns+1:2*ns, 2*ns+1:3*ns)     = res2_r * Phi;
    jac(2*ns+1:3*ns, 1:ns)          = res3_c_w * Phi;
    jac(2*ns+1:3*ns, ns+1:2*ns)     = res3_c_b * Phi;
    jac(2*ns+1:3*ns, 2*ns+1:3*ns)   = res3_r * Phi;
end

