function [res, jac, c_w, c_b, L, Y, I, q, Kp, Bp, r] = eval_resid_backward(c, c_0, param, glob, options, s, binding, mask)
%% Globals 
if nargin < 6
    s               = glob.s;
    s1              = glob.s;
    mask            = true(glob.Ns, 1);
else                                                                        % Solving period problem
    s1              = glob.s;
    chi_lev         = param.chi_lev;
end
if nargin < 7
    binding         = 'N';
end

ns              = size(s, 1);
ns1             = size(s1, 1);
basiscast       = glob.basiscast;
if options.disaster == 'Y'
    zind            = s(:, 3) == min(glob.zgrid);
    qual_shock      = (glob.dis_qual * (zind == 1) + 1 * (zind ~= 1));
else
    qual_shock      = ones(ns, 1);
end
K               = s(:, 1) .* qual_shock;
B               = min(s(:, 1) .* s(:, 2) - glob.transf, glob.bound_pol * s(:, 1));
Z               = s(:, 3);
Phi_Z           = splibas(glob.zgrid0, 0, glob.spliorder(3), s(:, 3));

if options.disaster == 'Y'
    zind            = s1(:, 3) == min(glob.zgrid);
    qual_shock      = (glob.dis_qual * (zind == 1) + 1 * (zind ~= 1));
else
    qual_shock      = ones(ns1, 1);
end
K1              = s1(:, 1) .* qual_shock;
B1              = min(s1(:, 1) .* s1(:, 2) - glob.transf, glob.bound_pol * s1(:, 1));
Z1              = s1(:, 3);

% Unpack
c1              = c(1:glob.Ns);
c2              = c(glob.Ns + 1:2 * glob.Ns);
c3              = c(2 * glob.Ns + 1:end);
c_w1            = c1;
c_b1            = c2;
r1              = c3;

c_01            = c_0(1:ns);
c_02            = c_0(ns + 1:2 * ns);
c_03            = c_0(2 * ns + 1:end);
c_w             = c_01;
c_b             = c_02;
r               = c_03;


%% Solve equations 
L               = solve_L_from_MRS(Z, K, c_w, param, glob, options);
Y               = production(Z, K, L, param, glob, options);
I               = Y - c_b - c_w;
q               = 1 ./ cap_prod_prime(I ./ K, param, glob, options);
Kp              = K .* (cap_prod(I ./ K, param, glob, options) + 1 - param.delta);
mpk             = production_k(Z, K, L, param, glob, options);
Pi              = cap_prod(I ./ K, param, glob, options) .* q - I ./ K;
Bp              = r .* (B + c_b + q .* Kp - mpk .* K - (1 - param.delta) * q .* K - Pi .* K);

% tomorrow
L1              = solve_L_from_MRS(Z1, K1, c_w1, param, glob, options);
Y1              = production(Z1, K1, L1, param, glob, options);
I1              = Y1 - c_b1 - c_w1;
q1              = 1 ./ cap_prod_prime(I1 ./ K1, param, glob, options);
Pi1             = cap_prod(I1 ./ K1, param, glob, options) .* q1 - I1 ./ K1;
mpk1            = production_k(Z1, K1, L1, param, glob, options);
mu_w1           = utility_c(c_w1, L1, param, glob, options);
mu_b1           = utility_c(c_b1, 0, param, glob, options);
mu_bprod1       = mu_b1 .* (mpk1 + (1 - param.delta) * q1 + Pi1) .* qual_shock;

% Enforce bounds
% Kp              = min(max(Kp, glob. kmin), glob.kmax);
% Bp              = min(max(Bp, glob. bmin), glob.bmax);

%% Create basis matrices for next states
ratiop          = Bp ./ Kp;

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder(2), ratiop);           % Basis for all b'
Phi_KBZp        = dprod(Phi_Z, dprod(Phi_Bp, Phi_Kp));                 % Basis for all (k', b', z')

%% Compute next state objects
PhiEmu_bprod    = basiscast * mu_bprod1;                                     % Approximate on the same basis, take expectation and get coefficients again
PhiEmu_b        = basiscast * mu_b1;                     
PhiEmu_w        = basiscast * mu_w1;                     
Emu_bprodp      = Phi_KBZp * PhiEmu_bprod;                               % Evaluate at (k', b', z)
Emu_bp          = Phi_KBZp * PhiEmu_b;                                   % Evaluate at (k', b', z)
Emu_wp          = Phi_KBZp * PhiEmu_w;                                   % Evaluate at (k', b', z)


%% Compute residuals
res                 = zeros(3 * ns, 1);
res(1:ns)           = (c_w) - utility_c_inv(glob.beta_w * r .* Emu_wp, L, param, glob, options);                    % Euler equation for workers
if options.lev_cap == 'Y' && binding == 'Y'
    res(ns+1:2*ns)      = chi_lev * q .* (utility_c(c_b, zeros(ns, 1), param, glob, options) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp)) - ... % Banker's euler equations
                        (utility_c(c_b, zeros(ns, 1), param, glob, options) .* q - (glob.beta_b * Emu_bprodp));
% 	res(ns+1:2*ns)      = chi_lev * q .* (ones(ns, 1) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp ./ utility_c(c_b, zeros(ns, 1), param, glob, options))) - ... % Banker's euler equations
%                         (q - (glob.beta_b * Emu_bprodp) ./ utility_c(c_b, zeros(ns, 1), param, glob, options));
    res(2*ns+1:3*ns)    = Bp - chi_lev * q .* Kp;                                                               % Binding leverage constraint
else
    res(ns+1:2*ns)      = (c_b) - utility_c_inv(glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp, zeros(ns, 1), param, glob, options);   % Euler equation for bankers
    res(2*ns+1:3*ns)    = utility_c_inv(Emu_bprodp, zeros(ns, 1), param, glob, options) - ...
                        utility_c_inv(Emu_bp .* q .* r ./ (1 - glob.dep_tax), zeros(ns, 1), param, glob, options);                      % Arbitrage
end

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
    
    Pi_c_w          = cap_prod(I ./ K, param, glob, options) .* q_c_w + cap_prod_prime(I ./ K, param, glob, options) .* I_c_w ./ K .* q - I_c_w ./ K;
    Pi_c_b          = cap_prod(I ./ K, param, glob, options) .* q_c_b + cap_prod_prime(I ./ K, param, glob, options) .* I_c_b ./ K .* q - I_c_b ./ K;
    Pi_r            = cap_prod(I ./ K, param, glob, options) .* q_r + cap_prod_prime(I ./ K, param, glob, options) .* I_r ./ K .* q - I_r ./ K;
    
    Bp_c_w          = r .* (q_c_w .* Kp + q .* Kp_c_w - mpk_c_w .* K - (1 - param.delta) * q_c_w .* K - Pi_c_w .* K);
    Bp_c_b          = r .* (ones(ns, 1) + q_c_b .* Kp + q .* Kp_c_b - mpk_c_b .* K - (1 - param.delta) * q_c_b .* K - Pi_c_b .* K);
    Bp_r            = ones(ns, 1) .* (B + c_b + q .* Kp - mpk .* K - (1 - param.delta) * q .* K - Pi .* K) + ...
                        r .* (q_r .* Kp + q .* Kp_r - mpk_r .* K - (1 - param.delta) * q_r .* K - Pi_r .* K);
                    
    %-------------
    Bp_c_w0         = Bp_c_w;
    Bp_c_b0         = Bp_c_b;
    Bp_r0           = Bp_r;
    Bp_c_w          = Bp_c_w ./ Kp - Bp ./ (Kp .^2) .* Kp_c_w;
    Bp_c_b          = Bp_c_b ./ Kp - Bp ./ (Kp .^2) .* Kp_c_b;
    Bp_r            = Bp_r ./ Kp - Bp ./ (Kp .^2) .* Kp_r;
    %-------------
    
    % Create derivative matrices
    Phi_Kp_der      = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp, 1);
    Phi_Bp_der      = splibas(glob.bgrid0, 0, glob.spliorder(2), ratiop, 1);
    Phi_KBZp_Kder   = dprod(Phi_Z, dprod(Phi_Bp, Phi_Kp_der));      
    Phi_KBZp_Bder   = dprod(Phi_Z, dprod(Phi_Bp_der, Phi_Kp));  
    
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
   
    Emu_bprodp_c_w  = Kder_bprod * diagKp_c_w + Bder_bprod * diagBp_c_w;
    Emu_bprodp_c_b  = Kder_bprod * diagKp_c_b + Bder_bprod * diagBp_c_b;
    Emu_bprodp_r    = Kder_bprod * diagKp_r + Bder_bprod * diagBp_r;
                    
    Emu_bp_c_w      = Kder_b * diagKp_c_w + Bder_b * diagBp_c_w;
    Emu_bp_c_b      = Kder_b * diagKp_c_b + Bder_b * diagBp_c_b;
    Emu_bp_r        = Kder_b * diagKp_r + Bder_b * diagBp_r;
    
    Emu_wp_c_w      = Kder_w * diagKp_c_w + Bder_w * diagBp_c_w;
    Emu_wp_c_b      = Kder_w * diagKp_c_b + Bder_w * diagBp_c_b;
    Emu_wp_r        = Kder_w * diagKp_r + Bder_w * diagBp_r;
    

    % Calculate the blocks of the jacobian
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
                
    if options.lev_cap == 'Y' && binding == 'Y'
        res2_c_w        = spdiags(chi_lev * q_c_w .* (utility_c(c_b, zeros(ns, 1), param, glob, options) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp)) + ...
                        chi_lev * q .* ( - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp_c_w)) - ...
                        (utility_c(c_b, zeros(ns, 1), param, glob, options) .* q_c_w - (glob.beta_b * Emu_bprodp_c_w)), 0, ns, ns);
        res2_c_b        = spdiags(chi_lev * q_c_b .* (utility_c(c_b, zeros(ns, 1), param, glob, options) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp)) + ...
                        chi_lev * q .* (utility_cc(c_b, zeros(ns, 1), param, glob, options) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp_c_b)) - ...
                        (utility_c(c_b, zeros(ns, 1), param, glob, options) .* q_c_b + ...
                        utility_cc(c_b, zeros(ns, 1), param, glob, options) .* q - (glob.beta_b * Emu_bprodp_c_b)), 0, ns, ns);
        res2_r          = spdiags(chi_lev * q_r .* (utility_c(c_b, zeros(ns, 1), param, glob, options) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp)) + ...
                        chi_lev * q .* ( - (glob.beta_b ./ (1 - glob.dep_tax) .* Emu_bp) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp_r)) - ...
                        r .* (utility_c(c_b, zeros(ns, 1), param, glob, options) .* q_r - (glob.beta_b * Emu_bprodp_r)), 0, ns, ns);
        
%         res2_c_w        = spdiags(chi_lev * q_c_w .* (ones(ns, 1) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp ./ utility_c(c_b, zeros(ns, 1), param, glob, options))) + ...
%                         chi_lev * q .* ( - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp_c_w ./ utility_c(c_b, zeros(ns, 1), param, glob, options))) - ...
%                         (q_c_w - (glob.beta_b * Emu_bprodp_c_w) ./ utility_c(c_b, zeros(ns, 1), param, glob, options)), 0, ns, ns);
%         
%         res2_c_b        = spdiags(chi_lev * q_c_b .* (ones(ns, 1) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp ./ utility_c(c_b, zeros(ns, 1), param, glob, options))) + ...
%                         chi_lev * q .* ( - glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp_c_b ./ utility_c(c_b, zeros(ns, 1), param, glob, options) + ... 
%                         (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp .* utility_cc(c_b, zeros(ns, 1), param, glob, options) ./ utility_c(c_b, zeros(ns, 1), param, glob, options) .^ 2))- ...
%                         (q_c_b - (glob.beta_b * Emu_bprodp_c_b) ./ utility_c(c_b, zeros(ns, 1), param, glob, options) + ...
%                         (glob.beta_b * Emu_bprodp .* utility_cc(c_b, zeros(ns, 1), param, glob, options)) ./ utility_c(c_b, zeros(ns, 1), param, glob, options) .^ 2), 0, ns, ns);
%                     
%         res2_r          = spdiags(chi_lev * q_r .* (ones(ns, 1) - (glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp ./ utility_c(c_b, zeros(ns, 1), param, glob, options))) + ...
%                         chi_lev * q .* ( - (glob.beta_b ./ (1 - glob.dep_tax) .* Emu_bp + glob.beta_b * r ./ (1 - glob.dep_tax) .* Emu_bp_r) ./ utility_c(c_b, zeros(ns, 1), param, glob, options)) - ...
%                         (q_r - (glob.beta_b * Emu_bprodp_r) ./ utility_c(c_b, zeros(ns, 1), param, glob, options)), 0, ns, ns);

        res3_c_w        = spdiags(Bp_c_w0 - chi_lev * q_c_w .* Kp - chi_lev * q .* Kp_c_w, 0, ns, ns);
        res3_c_b        = spdiags(Bp_c_b0 - chi_lev * q_c_b .* Kp - chi_lev * q .* Kp_c_b, 0, ns, ns);
        res3_r          = spdiags(Bp_r0 - chi_lev * q_r .* Kp - chi_lev * q .* Kp_r, 0, ns, ns);
    else
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
    end
    
    % Put together the jacobian
    jac                             = zeros(3 * ns, 3 * ns);
    jac(1:ns, 1:ns)                 = res1_c_w;
    jac(1:ns, ns+1:2*ns)            = res1_c_b;
    jac(1:ns, 2*ns+1:3*ns)          = res1_r;
    jac(ns+1:2*ns, 1:ns)            = res2_c_w;
    jac(ns+1:2*ns, ns+1:2*ns)       = res2_c_b;
    jac(ns+1:2*ns, 2*ns+1:3*ns)     = res2_r;
    jac(2*ns+1:3*ns, 1:ns)          = res3_c_w;
    jac(2*ns+1:3*ns, ns+1:2*ns)     = res3_c_b;
    jac(2*ns+1:3*ns, 2*ns+1:3*ns)   = res3_r;
end

