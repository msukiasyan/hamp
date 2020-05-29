function [res, jac, H_c_L, H_l_L, c_w, c_b, L, Y, I, q, Kp, Bp, gam] = eval_resid_FR(c, param, glob, options)
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
mu              = glob.s(:, 3);
Z               = glob.s(:, 4);
% Unpack
c1              = c(1:ns);
c2              = c(ns + 1:2 * ns);
c3              = c(2 * ns + 1:3 * ns);
% c4              = c(3 * ns + 1:4 * ns);
c_w1            = Phi * c1;
L1              = Phi * c2;
Bp1             = Phi * c3;
% c_b1            = Phi * c4;

c_01            = c(1:ns);
c_02            = c(ns + 1:2 * ns);
c_03            = c(2 * ns + 1:3 * ns);
% c_04            = c(3 * ns + 1:4 * ns);
c_w             = Phi * c_01;
L               = Phi * c_02;
Bp              = Phi * c_03;
% c_b             = Phi * c_04;


%% Solve equations 
H_l             = (B - c_w) .* utility_cl(c_w, L, param, glob, options) -  ...
    utility_ll(c_w, L, param, glob, options) .* L - utility_l(c_w, L, param, glob, options);
H_c             = (B - c_w) .* utility_cc(c_w, L, param, glob, options) -  ...
    utility_c(c_w, L, param, glob, options) - utility_cl(c_w, L, param, glob, options) .* L;

gam             = -((production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_l(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cl(c_w, L, param, glob, options)) ./ (-H_l - production_l(Z, K, L, param, glob, options) .* H_c));
c_b             = utility_c_inv(((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + gam .* H_c) / glob.lambda, 0, param, glob, options);
Y               = production(Z, K, L, param, glob, options);
I               = Y - c_b - c_w;
q               = 1 ./ cap_prod_prime(I ./ K, param, glob, options);
Kp              = K .* (cap_prod(I ./ K, param, glob, options) + 1 - param.delta);
mpk             = production_k(Z, K, L, param, glob, options);
mup             = gam;

% tomorrow
H1_l            = (B - c_w1) .* utility_cl(c_w1, L1, param, glob, options) -  ...
    utility_ll(c_w1, L1, param, glob, options) .* L1 - utility_l(c_w1, L1, param, glob, options);
H1_c            = (B - c_w1) .* utility_cc(c_w1, L1, param, glob, options) -  ...
    utility_c(c_w1, L1, param, glob, options) - utility_cl(c_w1, L1, param, glob, options) .* L1;

gam1            = -((production_l(Z, K, L1, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w1, L1, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w1, L1, param, glob, options)) + (1 - glob.lambda) * utility_l(c_w1, L1, param, glob, options) - ...
    (-mu) .* B .* utility_cl(c_w1, L1, param, glob, options)) ./ (-H1_l - production_l(Z, K, L1, param, glob, options) .* H1_c));
c_b1            = utility_c_inv(((1 - glob.lambda) * utility_c(c_w1, L1, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w1, L1, param, glob, options) + gam1 .* H1_c) / glob.lambda, 0, param, glob, options);

% Y1              = production(Z, K, L1, param, glob, options);
% I1              = Y1 - c_b1 - c_w1;
% q1              = 1 ./ cap_prod_prime(I1 ./ K, param, glob, options);
% mpk1            = production_k(Z, K, L1, param, glob, options);
Pi              = cap_prod(I ./ K, param, glob, options) .* q - I ./ K;
mu_w            = utility_c(c_w, L, param, glob, options);
mu_b            = utility_c(c_b, 0, param, glob, options);
mu_bprod        = mu_b .* (mpk + (1 - param.delta) * q + Pi);
mu_wgam         = mu_w .* gam;

% Enforce bounds
% Kp              = min(max(Kp, glob. kmin), glob.kmax);
% Bp              = min(max(Bp, glob. bmin), glob.bmax);

%% Create basis matrices for next states
ratiop          = Bp ./ Kp;

Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp);           % Basis for all k'
Phi_Bp          = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop);           % Basis for all b'
Phi_mup         = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup);           % Basis for all mu'
Phi_KBmuZp      = dprod(glob.Phi_Z, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp)));                 % Basis for all (k', b', z')

%% Compute next state objects
PhiEmu_bprod    = basiscast * mu_bprod;                                     % Approximate on the same basis, take expectation and get coefficients again
PhiEmu_wgam     = basiscast * mu_wgam;                     
PhiEmu_w        = basiscast * mu_w;                     
Emu_bprodp      = Phi_KBmuZp * PhiEmu_bprod;                               % Evaluate at (k', b', z)
Emu_wgamp       = Phi_KBmuZp * PhiEmu_wgam;                                   % Evaluate at (k', b', z)
Emu_wp          = Phi_KBmuZp * PhiEmu_w;                                   % Evaluate at (k', b', z)


%% Compute residuals
res                 = zeros(3 * ns, 1);
% res(1:ns)           = (B - c_w) .* utility_c(c_w, L, param, glob, options) - ...                                % LOM for deposits
%     utility_l(c_w, L, param, glob, options) .* L - glob.beta * Bp .* Emu_wp;                     
res(1:ns)           = (B - c_w) - ...                                % LOM for deposits
    utility_l(c_w, L, param, glob, options) .* L ./ utility_c(c_w, L, param, glob, options) ...
    - glob.beta * Bp .* Emu_wp ./ utility_c(c_w, L, param, glob, options);
res(ns+1:2*ns)      = (c_b) - utility_c_inv(glob.beta * Emu_bprodp ./ q, zeros(ns, 1), param, glob, options);   % Euler equation for bankers
res(2*ns+1:3*ns)    = (Emu_wgamp) - (gam .* Emu_wp);                                                                % Twisted martingale
% res(3*ns+1:4*ns)    = (c_b) - utility_c_inv(((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
%     (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + gam .* H_c) / glob.lambda, 0, param, glob, options); 

%% Compute the jacobian if requested
jac             = [];
if (nargout >= 2)
    % Compute derivatives with respect to all 4 inputs
    H_l_c_w         = (-1) .* utility_cl(c_w, L, param, glob, options) + ...
        (B - c_w) .* utility_ccl(c_w, L, param, glob, options) -  ...
        utility_cll(c_w, L, param, glob, options) .* L - utility_cl(c_w, L, param, glob, options);
%     H_l_c_b         = zeros(ns, 1);
    H_l_L           = (B - c_w) .* utility_cll(c_w, L, param, glob, options) -  ...
        utility_lll(c_w, L, param, glob, options) .* L - ...
        utility_ll(c_w, L, param, glob, options) - utility_ll(c_w, L, param, glob, options);
    H_l_Bp          = zeros(ns, 1);
    
    
    H_c_c_w         = (-1) .* utility_cc(c_w, L, param, glob, options) + ...
        (B - c_w) .* utility_ccc(c_w, L, param, glob, options) - ...
        utility_cc(c_w, L, param, glob, options) - utility_ccl(c_w, L, param, glob, options) .* L;
%     H_c_c_b         = zeros(ns, 1);
    H_c_L           = (B - c_w) .* utility_ccl(c_w, L, param, glob, options) -  ...
        utility_cl(c_w, L, param, glob, options) - (utility_cll(c_w, L, param, glob, options) .* L + ...
        utility_cl(c_w, L, param, glob, options));
    H_c_Bp          = zeros(ns, 1);
    
    gam_c_w         = -((production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_cc(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccc(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccl(c_w, L, param, glob, options)) .* (-H_l - production_l(Z, K, L, param, glob, options) .* H_c) + ...
        (-(-H_l_c_w - production_l(Z, K, L, param, glob, options) .* H_c_c_w)) .* (production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cc(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_l(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cl(c_w, L, param, glob, options))) ./ ((-H_l - production_l(Z, K, L, param, glob, options) .* H_c) .^ 2);
%     gam_c_b         = zeros(ns, 1);
    gam_L           = -((production_ll(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cc(c_w, L, param, glob, options)) + ...
        production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccl(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_ll(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cll(c_w, L, param, glob, options)) .* (-H_l - production_l(Z, K, L, param, glob, options) .* H_c) + ...
        (-(-H_l_L - production_l(Z, K, L, param, glob, options) .* H_c_L - production_ll(Z, K, L, param, glob, options) .* H_c)) .* (production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cc(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_l(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cl(c_w, L, param, glob, options))) ./ ((-H_l - production_l(Z, K, L, param, glob, options) .* H_c) .^ 2);
    gam_Bp          = zeros(ns, 1);
    
    c_b_c_w         = utility_c_inv_prime(((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + gam .* H_c) / glob.lambda, 0, param, glob, options) .* ...
        (((1 - glob.lambda) * utility_cc(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccc(c_w, L, param, glob, options) + gam_c_w .* H_c + gam .* H_c_c_w) / glob.lambda);
    c_b_L           = utility_c_inv_prime(((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + gam .* H_c) / glob.lambda, 0, param, glob, options) .* ...
        (((1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccl(c_w, L, param, glob, options) + gam_L .* H_c + gam .* H_c_L) / glob.lambda);
    c_b_Bp          = zeros(ns, 1);
    
    Y_c_w           = zeros(ns, 1);
%     Y_c_b           = zeros(ns, 1);
    Y_L             = production_l(Z, K, L, param, glob, options);
    Y_Bp            = zeros(ns, 1);
    
    I_c_w           = Y_c_w - ones(ns, 1) - c_b_c_w;
%     I_c_b           = Y_c_b - ones(ns, 1);
    I_L             = Y_L - c_b_L;
    I_Bp            = Y_Bp - c_b_Bp;
    
    q_c_w           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c_w;
%     q_c_b           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c_b;
    q_L             = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_L;
    q_Bp            = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_Bp;
    
    Kp_c_w          = (1 ./ q) .* I_c_w;
%     Kp_c_b          = (1 ./ q) .* I_c_b;
    Kp_L            = (1 ./ q) .* I_L;
    Kp_Bp           = (1 ./ q) .* I_Bp;
                    
    %-------------
    Bp_c_w          = - Bp ./ (Kp .^2) .* Kp_c_w;
%     Bp_c_b          = - Bp ./ (Kp .^2) .* Kp_c_b;
    Bp_L            = - Bp ./ (Kp .^2) .* Kp_L;
    Bp_Bp           = 1 ./ Kp - Bp ./ (Kp .^2) .* Kp_Bp;
    %-------------
    
    % Create derivative matrices
    Phi_Kp_der      = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp, 1);
    Phi_Bp_der      = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop, 1);
    Phi_mup_der     = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup, 1);
    Phi_KBmuZp_Kder     = dprod(glob.Phi_Z, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp_der)));
    Phi_KBmuZp_Bder     = dprod(glob.Phi_Z, dprod(Phi_mup, dprod(Phi_Bp_der, Phi_Kp)));
    Phi_KBmuZp_muder    = dprod(glob.Phi_Z, dprod(Phi_mup_der, dprod(Phi_Bp, Phi_Kp)));
    
    Kder_bprod      = spdiags(Phi_KBmuZp_Kder * PhiEmu_bprod, 0, ns, ns);
    Bder_bprod      = spdiags(Phi_KBmuZp_Bder * PhiEmu_bprod, 0, ns, ns);
    muder_bprod     = spdiags(Phi_KBmuZp_muder * PhiEmu_bprod, 0, ns, ns);
    Kder_w          = spdiags(Phi_KBmuZp_Kder * PhiEmu_w, 0, ns, ns);
    Bder_w          = spdiags(Phi_KBmuZp_Bder * PhiEmu_w, 0, ns, ns);
    muder_w         = spdiags(Phi_KBmuZp_muder * PhiEmu_w, 0, ns, ns);
    Kder_wgam       = spdiags(Phi_KBmuZp_Kder * PhiEmu_wgam, 0, ns, ns);
    Bder_wgam       = spdiags(Phi_KBmuZp_Bder * PhiEmu_wgam, 0, ns, ns);
    muder_wgam      = spdiags(Phi_KBmuZp_muder * PhiEmu_wgam, 0, ns, ns);
    
    diagKp_c_w      = spdiags(Kp_c_w, 0, ns, ns);
%     diagKp_c_b      = spdiags(Kp_c_b, 0, ns, ns);
    diagKp_L        = spdiags(Kp_L, 0, ns, ns);
    diagKp_Bp       = spdiags(Kp_Bp, 0, ns, ns);
    diagBp_c_w      = spdiags(Bp_c_w, 0, ns, ns);
%     diagBp_c_b      = spdiags(Bp_c_b, 0, ns, ns);
    diagBp_L        = spdiags(Bp_L, 0, ns, ns);
    diagBp_Bp       = spdiags(Bp_Bp, 0, ns, ns);
    diagmup_c_w     = spdiags(gam_c_w, 0, ns, ns);
%     diagmup_c_b     = spdiags(gam_c_b, 0, ns, ns);
    diagmup_L       = spdiags(gam_L, 0, ns, ns);
    diagmup_Bp      = spdiags(gam_Bp, 0, ns, ns);
    
    
    next_basis      = Phi_KBmuZp * basiscast;
    
    mu_w_c_w        = utility_cc(c_w, L, param, glob, options) .* ones(ns, 1);
%     mu_w_c_b        = zeros(ns, 1);
    mu_w_L          = utility_cl(c_w, L, param, glob, options);
    mu_w_Bp         = zeros(ns, 1);
    
    mu_b_c_w        = utility_cc(c_b, 0, param, glob, options) .* c_b_c_w;
%     mu_b_c_b        = utility_cc(c_b, 0, param, glob, options) .* ones(ns, 1);
    mu_b_L          = utility_cc(c_b, 0, param, glob, options) .* c_b_L;
    mu_b_Bp         = utility_cc(c_b, 0, param, glob, options) .* c_b_Bp;
    
    mpk_c_w         = zeros(ns, 1);
%     mpk_c_b         = zeros(ns, 1);
    mpk_L           = production_kl(Z, K, L, param, glob, options);
    mpk_Bp          = zeros(ns, 1);
    
    Pi_c_w          = cap_prod_prime(I ./ K, param, glob, options) .* (I_c_w ./ K) .* q + cap_prod(I ./ K, param, glob, options) .* q_c_w - ...
        I_c_w ./ K;
%     Pi_c_b          = cap_prod_prime(I ./ K, param, glob, options) .* (I_c_b ./ K) .* q + cap_prod(I ./ K, param, glob, options) .* q_c_b - ...
%         I_c_b ./ K;
    Pi_L            = cap_prod_prime(I ./ K, param, glob, options) .* (I_L ./ K) .* q + cap_prod(I ./ K, param, glob, options) .* q_L - ...
        I_L ./ K;
    Pi_Bp           = cap_prod_prime(I ./ K, param, glob, options) .* (I_Bp ./ K) .* q + cap_prod(I ./ K, param, glob, options) .* q_Bp - ...
        I_Bp ./ K;
    
    mu_bprod_c_w    = mu_b_c_w .* (mpk + (1 - param.delta) * q + Pi) + mu_b .* (mpk_c_w + (1 - param.delta) * q_c_w + Pi_c_w);
%     mu_bprod_c_b    = mu_b_c_b .* (mpk + (1 - param.delta) * q + Pi) + mu_b .* (mpk_c_b + (1 - param.delta) * q_c_b + Pi_c_b);
    mu_bprod_L      = mu_b_L .* (mpk + (1 - param.delta) * q + Pi) + mu_b .* (mpk_L + (1 - param.delta) * q_L + Pi_L);
    mu_bprod_Bp     = mu_b_Bp .* (mpk + (1 - param.delta) * q + Pi) + mu_b .* (mpk_Bp + (1 - param.delta) * q_Bp + Pi_Bp);
    
    mu_wgam_c_w     = mu_w_c_w .* gam + mu_w .* gam_c_w;
%     mu_wgam_c_b     = mu_w_c_b .* gam + mu_w .* gam_c_b;
    mu_wgam_L       = mu_w_L .* gam + mu_w .* gam_L;
    mu_wgam_Bp      = mu_w_Bp .* gam + mu_w .* gam_Bp;
    
   
    Emu_bprodp_c_w  = Kder_bprod * diagKp_c_w + Bder_bprod * diagBp_c_w + muder_bprod * diagmup_c_w + ...
        next_basis * spdiags(mu_bprod_c_w, 0, ns, ns);
%     Emu_bprodp_c_b  = Kder_bprod * diagKp_c_b + Bder_bprod * diagBp_c_b + muder_bprod * diagmup_c_b + ...
%         next_basis * spdiags(mu_bprod_c_b, 0, ns, ns);
    Emu_bprodp_L    = Kder_bprod * diagKp_L + Bder_bprod * diagBp_L + muder_bprod * diagmup_L + ...
        next_basis * spdiags(mu_bprod_L, 0, ns, ns);
    Emu_bprodp_Bp   = Kder_bprod * diagKp_Bp + Bder_bprod * diagBp_Bp + muder_bprod * diagmup_Bp + ...
        next_basis * spdiags(mu_bprod_Bp, 0, ns, ns);
    
    Emu_wp_c_w      = Kder_w * diagKp_c_w + Bder_w * diagBp_c_w + muder_w * diagmup_c_w + ...
        next_basis * spdiags(mu_w_c_w, 0, ns, ns);
%     Emu_wp_c_b      = Kder_w * diagKp_c_b + Bder_w * diagBp_c_b + muder_w * diagmup_c_b + ...
%         next_basis * spdiags(mu_w_c_b, 0, ns, ns);
    Emu_wp_L        = Kder_w * diagKp_L + Bder_w * diagBp_L + muder_w * diagmup_L + ...
        next_basis * spdiags(mu_w_L, 0, ns, ns);
    Emu_wp_Bp       = Kder_w * diagKp_Bp + Bder_w * diagBp_Bp + muder_w * diagmup_Bp + ...
        next_basis * spdiags(mu_w_Bp, 0, ns, ns);
                    
    Emu_wgamp_c_w   = Kder_wgam * diagKp_c_w + Bder_wgam * diagBp_c_w + muder_wgam * diagmup_c_w + ...
        next_basis * spdiags(mu_wgam_c_w, 0, ns, ns);
%     Emu_wgamp_c_b   = Kder_wgam * diagKp_c_b + Bder_wgam * diagBp_c_b + muder_wgam * diagmup_c_b + ...
%         next_basis * spdiags(mu_wgam_c_b, 0, ns, ns);
    Emu_wgamp_L     = Kder_wgam * diagKp_L + Bder_wgam * diagBp_L + muder_wgam * diagmup_L + ...
        next_basis * spdiags(mu_wgam_L, 0, ns, ns);
    Emu_wgamp_Bp    = Kder_wgam * diagKp_Bp + Bder_wgam * diagBp_Bp + muder_wgam * diagmup_Bp + ...
        next_basis * spdiags(mu_wgam_Bp, 0, ns, ns);
    

    % Calculate the blocks of the jacobian
    % mu_prime1       = spdiags(utility_c_inv_prime(glob.beta_w * r .* Emu_wp, L, param, glob, options), 0, ns, ns);
    mu_prime2       = spdiags(utility_c_inv_prime(glob.beta * Emu_bprodp ./ q, zeros(ns, 1), param, glob, options), 0, ns, ns);
    mu_prime4       = spdiags(utility_c_inv_prime(((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
            (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + gam .* H_c) / glob.lambda, 0, param, glob, options), 0, ns, ns);
    mu_prime31      = spdiags(utility_c_inv_prime(Emu_wgamp, 0, param, glob, options), 0, ns, ns);
    mu_prime32      = spdiags(utility_c_inv_prime(gam .* Emu_wp, 0, param, glob, options), 0, ns, ns);
    
%     res1_c_w        = spdiags((-1) .* utility_c(c_w, L, param, glob, options) + ...
%         (B - c_w) .* utility_cc(c_w, L, param, glob, options) - ...
%         utility_cl(c_w, L, param, glob, options) .* L, 0, ns, ns) - glob.beta * spdiags(Bp, 0, ns, ns) * Emu_wp_c_w;
%     res1_c_b        = - glob.beta * spdiags(Bp, 0, ns, ns) * Emu_wp_c_b;
%     res1_L          = spdiags((B - c_w) .* utility_cl(c_w, L, param, glob, options) - ...     
%         utility_ll(c_w, L, param, glob, options) .* L - utility_l(c_w, L, param, glob, options), 0, ns, ns) - ...
%         glob.beta * spdiags(Bp, 0, ns, ns) * Emu_wp_L;
%     res1_Bp         = - glob.beta * (spdiags(Bp, 0, ns, ns) * Emu_wp_Bp + spdiags(Emu_wp, 0, ns, ns) * speye(ns, ns));

    res(1:ns)           = (B - c_w) - ...                                % LOM for deposits
        utility_l(c_w, L, param, glob, options) .* L ./ utility_c(c_w, L, param, glob, options) ...
        - glob.beta * Bp .* Emu_wp ./ utility_c(c_w, L, param, glob, options);

    res1_c_w        = spdiags((-1) - utility_cl(c_w, L, param, glob, options) .* L ./ utility_c(c_w, L, param, glob, options) + ...
        utility_cc(c_w, L, param, glob, options) .* utility_l(c_w, L, param, glob, options) .* L ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns) - ...
        glob.beta * spdiags(Bp ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * Emu_wp_c_w + ...
        glob.beta * spdiags(utility_cc(c_w, L, param, glob, options) .* Bp .* Emu_wp ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns);
%     res1_c_b        = - glob.beta * spdiags(Bp ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * Emu_wp_c_b;
    res1_L          = spdiags((B - c_w) - ...
        (utility_ll(c_w, L, param, glob, options) .* L + utility_l(c_w, L, param, glob, options)) ./ utility_c(c_w, L, param, glob, options) + ...
        utility_cl(c_w, L, param, glob, options) .* utility_l(c_w, L, param, glob, options) .* L ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns) - ...
        glob.beta * spdiags(Bp ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * Emu_wp_L + ...
        glob.beta * spdiags(Bp .* Emu_wp .* utility_cl(c_w, L, param, glob, options) ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns);
    res1_Bp         = - glob.beta * (spdiags(Bp ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * Emu_wp_Bp + ... 
        spdiags(Emu_wp ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * speye(ns, ns));
    
    res2_c_w        = spdiags(c_b_c_w, 0, ns, ns) - mu_prime2 * glob.beta * ...
        (Emu_bprodp_c_w * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_c_w .* Emu_bprodp ./ (q .^ 2), 0, ns, ns));
%     res2_c_b        = speye(ns, ns) - mu_prime2 * ...
%         (Emu_bprodp_c_b * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_c_b .* Emu_bprodp ./ (q .^ 2), 0, ns, ns));
    res2_L          = spdiags(c_b_L, 0, ns, ns) - mu_prime2 * glob.beta * ...
        (Emu_bprodp_L * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_L .* Emu_bprodp ./ (q .^ 2), 0, ns, ns));
    res2_Bp         = spdiags(c_b_Bp, 0, ns, ns) - mu_prime2 * glob.beta * ...
        (Emu_bprodp_Bp * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_Bp .* Emu_bprodp ./ (q .^ 2), 0, ns, ns));
    
    res3_c_w        = Emu_wgamp_c_w - (spdiags(gam, 0, ns, ns) * Emu_wp_c_w + spdiags(gam_c_w .* Emu_wp, 0, ns, ns));
%     res3_c_b        = Emu_wgamp_c_b - (spdiags(gam, 0, ns, ns) * Emu_wp_c_b + spdiags(gam_c_b .* Emu_wp, 0, ns, ns));
    res3_L          = Emu_wgamp_L - (spdiags(gam, 0, ns, ns) * Emu_wp_L + spdiags(gam_L .* Emu_wp, 0, ns, ns));
    res3_Bp         = Emu_wgamp_Bp - (spdiags(gam, 0, ns, ns) * Emu_wp_Bp + spdiags(gam_Bp .* Emu_wp, 0, ns, ns));
    
%     res4_c_w        = mu_prime4 * spdiags(- ( ((1 - glob.lambda) * utility_cc(c_w, L, param, glob, options) - ...
%         (-mu) .* B .* utility_ccc(c_w, L, param, glob, options) + gam .* H_c_c_w + gam_c_w .* H_c) / glob.lambda ), 0, ns, ns);
%     % res4_c_b        = spdiags(utility_cc(c_b, 0, param, glob, options) - ( (gam .* H_c_c_b + gam_c_b .* H_c) / glob.lambda ), 0, ns, ns);
%     res4_c_b        = speye(ns, ns) - mu_prime4 * spdiags(( (gam .* H_c_c_b + gam_c_b .* H_c) / glob.lambda ), 0, ns, ns);
%     res4_L          = mu_prime4 * spdiags(- (((1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - ...
%         (-mu) .* B .* utility_ccl(c_w, L, param, glob, options) + gam_L .* H_c + gam .* H_c_L) / glob.lambda), 0, ns, ns);
%     res4_Bp         = mu_prime4 * spdiags(- ((gam .* H_c_Bp + gam_Bp .* H_c) / glob.lambda), 0, ns, ns);
    
    % Put together the jacobian
%     jac                             = (zeros(4 * ns, 4 * ns));
%     jac(1:ns, 1:ns)                 = res1_c_w * Phi;
%     jac(1:ns, ns+1:2*ns)            = res1_L * Phi;
%     jac(1:ns, 2*ns+1:3*ns)          = res1_Bp * Phi;
%     jac(1:ns, 3*ns+1:4*ns)          = res1_c_b * Phi;
%     jac(ns+1:2*ns, 1:ns)            = res2_c_w * Phi;
%     jac(ns+1:2*ns, ns+1:2*ns)       = res2_L * Phi;
%     jac(ns+1:2*ns, 2*ns+1:3*ns)     = res2_Bp * Phi;
%     jac(ns+1:2*ns, 3*ns+1:4*ns)     = res2_c_b * Phi;
%     jac(2*ns+1:3*ns, 1:ns)          = res3_c_w * Phi;
%     jac(2*ns+1:3*ns, ns+1:2*ns)     = res3_L * Phi;
%     jac(2*ns+1:3*ns, 2*ns+1:3*ns)   = res3_Bp * Phi;
%     jac(2*ns+1:3*ns, 3*ns+1:4*ns)   = res3_c_b * Phi;
%     jac(3*ns+1:4*ns, 1:ns)          = res4_c_w * Phi;
%     jac(3*ns+1:4*ns, ns+1:2*ns)     = res4_L * Phi;
%     jac(3*ns+1:4*ns, 2*ns+1:3*ns)   = res4_Bp * Phi;
%     jac(3*ns+1:4*ns, 3*ns+1:4*ns)   = res4_c_b * Phi;
    jac                             = (zeros(3 * ns, 3 * ns));
    jac(1:ns, 1:ns)                 = res1_c_w * Phi;
    jac(1:ns, ns+1:2*ns)            = res1_L * Phi;
    jac(1:ns, 2*ns+1:3*ns)          = res1_Bp * Phi;
    jac(ns+1:2*ns, 1:ns)            = res2_c_w * Phi;
    jac(ns+1:2*ns, ns+1:2*ns)       = res2_L * Phi;
    jac(ns+1:2*ns, 2*ns+1:3*ns)     = res2_Bp * Phi;
    jac(2*ns+1:3*ns, 1:ns)          = res3_c_w * Phi;
    jac(2*ns+1:3*ns, ns+1:2*ns)     = res3_L * Phi;
    jac(2*ns+1:3*ns, 2*ns+1:3*ns)   = res3_Bp * Phi;
end

