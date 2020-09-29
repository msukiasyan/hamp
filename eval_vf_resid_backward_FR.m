function [res, jac, eq] = eval_vf_resid_backward_FR(EWc, c_0, param, glob, options)
%% Unpack
s               = glob.s;  
ns              = size(s, 1);
Phip            = glob.Phisp;
Phi             = eye(ns, ns);
Phi_Z           = glob.Phi_Z;

if options.disaster == 'Y'
    zind            = s(:, 4) == min(glob.zgrid);
    qual_shock      = (glob.dis_qual * (zind == 1) + 1 * (zind ~= 1));
else
    qual_shock      = ones(ns, 1);
end
K               = glob.s(:, 1) .* qual_shock;
B               = glob.s(:, 1) .* glob.s(:, 2);
mu              = glob.s(:, 3);
Z               = glob.s(:, 4);

c_w             = c_0(1:ns);
L               = c_0(ns + 1:2 * ns);
Bp              = c_0(2 * ns + 1:3 * ns);
c_b             = c_0(3 * ns + 1:4 * ns);
gam             = c_0(4 * ns + 1:5 * ns);

if strcmp(options.FB, 'Y')
    gam         = zeros(ns, 1);
    mu          = zeros(ns, 1);
    Bp          = glob.s(:, 1) .* glob.s(:, 2);
end


%% Solve equations 
H_l             = (B - c_w) .* utility_cl(c_w, L, param, glob, options) -  ...
    utility_ll(c_w, L, param, glob, options) .* L - utility_l(c_w, L, param, glob, options);
H_c             = (B - c_w) .* utility_cc(c_w, L, param, glob, options) -  ...
    utility_c(c_w, L, param, glob, options) - utility_cl(c_w, L, param, glob, options) .* L;

Y               = production(Z, K, L, param, glob, options);
I               = Y - c_b - c_w;
q               = 1 ./ cap_prod_prime(I ./ K, param, glob, options);
Kp              = K .* (cap_prod(I ./ K, param, glob, options) + 1 - param.delta);
mpk             = production_k(Z, K, L, param, glob, options);
mup             = gam;
ratiop          = Bp ./ Kp;

%% Create derivative matrices
derPhi_K0       = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp, 0);  
derPhi_K1       = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp, 1);  
derPhi_rat0     = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop, 0);  
derPhi_rat1     = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop, 1);  
derPhi_mu0      = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup, 0);  
derPhi_mu1      = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup, 1);  

%% Next period derivatives
derK            = dprod(Phi_Z, dprod(derPhi_mu0, dprod(derPhi_rat0, derPhi_K1)));
derrat          = dprod(Phi_Z, dprod(derPhi_mu0, dprod(derPhi_rat1, derPhi_K0)));
dermu           = dprod(Phi_Z, dprod(derPhi_mu1, dprod(derPhi_rat0, derPhi_K0)));

%% Compute expected continuation values
EW_K            = derK * EWc;
EW_B            = derrat * EWc;
EW_mu           = dermu * EWc;

%% Compute residuals
res                 = zeros(3 * ns, 1);
% Constraint
res(1:ns)           = (B - c_w) - ...                                                                                   
    utility_l(c_w, L, param, glob, options) .* L ./ utility_c(c_w, L, param, glob, options) ...
    + glob.beta * (-EW_mu) ./ utility_c(c_w, L, param, glob, options);
% Intertemporal decision
res(ns+1:2*ns)      = (c_b) - utility_c_inv(glob.beta * EW_K ./ q / glob.lambda, zeros(ns, 1), param, glob, options);
% FOC wrt B
res(2*ns+1:3*ns)    = EW_B ./ K;
% Marginal utility equalisation
res(3*ns+1:4*ns)    = glob.lambda * utility_c(c_b, 0, param, glob, options) - ...
    ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - (-mu) .* B .* utility_cc(c_w, L, param, glob, options) + (-gam) .* H_c);
% Labor decision
res(4*ns+1:5*ns)    = (-gam) .* (-H_l - production_l(Z, K, L, param, glob, options) .* H_c) - (production_l(Z, K, L, param, glob, options) .* ((1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cc(c_w, L, param, glob, options)) + (1 - glob.lambda) * utility_l(c_w, L, param, glob, options) - ...
    (-mu) .* B .* utility_cl(c_w, L, param, glob, options));

if strcmp(options.FB, 'Y')
    res(1:ns)           = zeros(ns, 1);
    res(2*ns+1:3*ns)    = zeros(ns, 1);
end

%% Pack
eq.c_w                  = c_w;
eq.c_b                  = c_b;
eq.L                    = L;
eq.Kp                   = Kp;
eq.Bp                   = Bp;
eq.mup                  = mup;
eq.Y                    = Y;
eq.I                    = I;
eq.q                    = q;

%% Compute the jacobian if requested
jac             = [];
if (nargout == 2)
    % Compute derivatives with respect to all 5 inputs;
    % lots of chain rules...
    H_l_c_w         = (-1) .* utility_cl(c_w, L, param, glob, options) + ...
        (B - c_w) .* utility_ccl(c_w, L, param, glob, options) -  ...
        utility_cll(c_w, L, param, glob, options) .* L - utility_cl(c_w, L, param, glob, options);
    H_l_L           = (B - c_w) .* utility_cll(c_w, L, param, glob, options) -  ...
        utility_lll(c_w, L, param, glob, options) .* L - ...
        utility_ll(c_w, L, param, glob, options) - utility_ll(c_w, L, param, glob, options);
    H_l_Bp          = zeros(ns, 1);
    
    
    H_c_c_w         = (-1) .* utility_cc(c_w, L, param, glob, options) + ...
        (B - c_w) .* utility_ccc(c_w, L, param, glob, options) - ...
        utility_cc(c_w, L, param, glob, options) - utility_ccl(c_w, L, param, glob, options) .* L;
    H_c_L           = (B - c_w) .* utility_ccl(c_w, L, param, glob, options) -  ...
        utility_cl(c_w, L, param, glob, options) - (utility_cll(c_w, L, param, glob, options) .* L + ...
        utility_cl(c_w, L, param, glob, options));
    H_c_Bp          = zeros(ns, 1);

    gam_c_w         = zeros(ns, 1);
    gam_L           = zeros(ns, 1);
    gam_Bp          = zeros(ns, 1);
    gam_c_b         = zeros(ns, 1);
    gam_gam         = ones(ns, 1);
    
    c_b_c_w         = zeros(ns, 1);
    c_b_L           = zeros(ns, 1);
    c_b_Bp          = zeros(ns, 1);
    c_b_c_b         = ones(ns, 1);
    c_b_gam         = zeros(ns, 1);
    
    Y_c_w           = zeros(ns, 1);
    Y_L             = production_l(Z, K, L, param, glob, options);
    Y_Bp            = zeros(ns, 1);
    Y_gam           = zeros(ns, 1);
    
    I_c_w           = Y_c_w - ones(ns, 1) - c_b_c_w;
    I_L             = Y_L - c_b_L;
    I_Bp            = Y_Bp - c_b_Bp;
    I_c_b           = - c_b_c_b;
    I_gam           = zeros(ns, 1);
    
    q_c_w           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c_w;
    q_L             = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_L;
    q_Bp            = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_Bp;
    q_c_b           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c_b;
    q_gam           = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_gam;
    
    Kp_c_w          = (1 ./ q) .* I_c_w;
    Kp_L            = (1 ./ q) .* I_L;
    Kp_Bp           = (1 ./ q) .* I_Bp;
    Kp_c_b          = (1 ./ q) .* I_c_b;
    Kp_gam          = (1 ./ q) .* I_gam;
                    
    Bp_c_w          = - Bp ./ (Kp .^2) .* Kp_c_w;
    Bp_L            = - Bp ./ (Kp .^2) .* Kp_L;
    Bp_Bp           = 1 ./ Kp - Bp ./ (Kp .^2) .* Kp_Bp;
    Bp_c_b          = - Bp ./ (Kp .^2) .* Kp_c_b;
    Bp_gam          = - Bp ./ (Kp .^2) .* Kp_gam;
    
    if strcmp(options.FB, 'Y')
        Bp_Bp       = - Bp ./ (Kp .^2) .* Kp_Bp;
    end
    
    % Create second-order derivative matrices
    derPhi_K2       = splibas(glob.kgrid0, 0, glob.spliorder_FR(1), Kp, 2); 
    derPhi_rat2     = splibas(glob.bgrid0, 0, glob.spliorder_FR(2), ratiop, 2);
    derPhi_mu2      = splibas(glob.mugrid0, 0, glob.spliorder_FR(3), mup, 2);
    
    derKK           = spdiags(dprod(Phi_Z, dprod(derPhi_mu0, dprod(derPhi_rat0, derPhi_K2))) * EWc, 0, ns, ns);
    derKrat         = spdiags(dprod(Phi_Z, dprod(derPhi_mu0, dprod(derPhi_rat1, derPhi_K1))) * EWc, 0, ns, ns);
    derKmu          = spdiags(dprod(Phi_Z, dprod(derPhi_mu1, dprod(derPhi_rat0, derPhi_K1))) * EWc, 0, ns, ns);
    derratrat       = spdiags(dprod(Phi_Z, dprod(derPhi_mu0, dprod(derPhi_rat2, derPhi_K0))) * EWc, 0, ns, ns);
    derratmu        = spdiags(dprod(Phi_Z, dprod(derPhi_mu1, dprod(derPhi_rat1, derPhi_K0))) * EWc, 0, ns, ns);
    dermumu         = spdiags(dprod(Phi_Z, dprod(derPhi_mu2, dprod(derPhi_rat0, derPhi_K0))) * EWc, 0, ns, ns);
    
    diagKp_c_w      = spdiags(Kp_c_w, 0, ns, ns);
    diagKp_L        = spdiags(Kp_L, 0, ns, ns);
    diagKp_Bp       = spdiags(Kp_Bp, 0, ns, ns);
    diagKp_c_b      = spdiags(Kp_c_b, 0, ns, ns);
    
    diagBp_c_w      = spdiags(Bp_c_w, 0, ns, ns);
    diagBp_L        = spdiags(Bp_L, 0, ns, ns);
    diagBp_Bp       = spdiags(Bp_Bp, 0, ns, ns);
    diagBp_c_b      = spdiags(Bp_c_b, 0, ns, ns);
    
    diagmup_c_w     = spdiags(gam_c_w, 0, ns, ns);
    diagmup_L       = spdiags(gam_L, 0, ns, ns);
    diagmup_Bp      = spdiags(gam_Bp, 0, ns, ns);
    diagmup_c_b     = spdiags(gam_c_b, 0, ns, ns);
   
    EW_K_c_w        = derKK * diagKp_c_w + derKrat * diagBp_c_w + derKmu * diagmup_c_w;
    EW_K_L          = derKK * diagKp_L + derKrat * diagBp_L + derKmu * diagmup_L;
    EW_K_Bp         = derKK * diagKp_Bp + derKrat * diagBp_Bp + derKmu * diagmup_Bp;
    EW_K_c_b        = derKK * diagKp_c_b + derKrat * diagBp_c_b + derKmu * diagmup_c_b;
    EW_K_gam        = derKmu;
    
    EW_B_c_w        = derKrat * diagKp_c_w + derratrat * diagBp_c_w + derratmu * diagmup_c_w;
    EW_B_L          = derKrat * diagKp_L + derratrat * diagBp_L + derratmu * diagmup_L;
    EW_B_Bp         = derKrat * diagKp_Bp + derratrat  * diagBp_Bp + derratmu * diagmup_Bp;
    EW_B_c_b        = derKrat * diagKp_c_b + derratrat * diagBp_c_b + derratmu * diagmup_c_b;
    EW_B_gam        = derratmu;
    
    EW_mu_c_w       = derKmu * diagKp_c_w + derratmu * diagBp_c_w + dermumu * diagmup_c_w;
    EW_mu_L         = derKmu * diagKp_L + derratmu * diagBp_L + dermumu * diagmup_L;
    EW_mu_Bp        = derKmu * diagKp_Bp + derratmu * diagBp_Bp + dermumu * diagmup_Bp;
    EW_mu_c_b       = derKmu * diagKp_c_b + derratmu * diagBp_c_b + dermumu * diagmup_c_b;
    EW_mu_gam       = dermumu;
                    
    % Calculate the blocks of the jacobian
    mu_prime2       = spdiags(utility_c_inv_prime(glob.beta / glob.lambda * EW_K ./ q, zeros(ns, 1), param, glob, options), 0, ns, ns);
    
    res1_c_w        = spdiags((-1) - utility_cl(c_w, L, param, glob, options) .* L ./ utility_c(c_w, L, param, glob, options) + ...
        utility_cc(c_w, L, param, glob, options) .* utility_l(c_w, L, param, glob, options) .* L ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns) + ...
        glob.beta * spdiags(1 ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * (-EW_mu_c_w) - ...
        glob.beta * spdiags(utility_cc(c_w, L, param, glob, options) .* (-EW_mu) ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns);
    res1_L          = spdiags((B - c_w) - ...
        (utility_ll(c_w, L, param, glob, options) .* L + utility_l(c_w, L, param, glob, options)) ./ utility_c(c_w, L, param, glob, options) + ...
        utility_cl(c_w, L, param, glob, options) .* utility_l(c_w, L, param, glob, options) .* L ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns) + ...
        glob.beta * spdiags(1 ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * (-EW_mu_L) - ...
        glob.beta * spdiags((-EW_mu) .* utility_cl(c_w, L, param, glob, options) ./ (utility_c(c_w, L, param, glob, options) .^ 2), 0, ns, ns);
    res1_Bp         = glob.beta * spdiags(1 ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * (-EW_mu_Bp);
    res1_c_b        = glob.beta * spdiags(1 ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * (-EW_mu_c_b);
    res1_gam        = glob.beta * spdiags(1 ./ utility_c(c_w, L, param, glob, options), 0, ns, ns) * (-EW_mu_gam);
    
    res2_c_w        = spdiags(c_b_c_w, 0, ns, ns) - mu_prime2 * glob.beta / glob.lambda * ...
        (EW_K_c_w * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_c_w .* EW_K ./ (q .^ 2), 0, ns, ns));
    res2_L          = spdiags(c_b_L, 0, ns, ns) - mu_prime2 * glob.beta / glob.lambda * ...
        (EW_K_L * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_L .* EW_K ./ (q .^ 2), 0, ns, ns));
    res2_Bp         = spdiags(c_b_Bp, 0, ns, ns) - mu_prime2 * glob.beta / glob.lambda * ...
        (EW_K_Bp * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_Bp .* EW_K ./ (q .^ 2), 0, ns, ns));
    res2_c_b        = spdiags(c_b_c_b, 0, ns, ns) - mu_prime2 * glob.beta / glob.lambda * ...
        (EW_K_c_b * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_c_b .* EW_K ./ (q .^ 2), 0, ns, ns));
    res2_gam        = spdiags(c_b_gam, 0, ns, ns) - mu_prime2 * glob.beta / glob.lambda * ...
        (EW_K_gam * spdiags(1 ./ q, 0, ns, ns) + spdiags(-q_gam .* EW_K ./ (q .^ 2), 0, ns, ns));
    
    res3_c_w        = EW_B_c_w * spdiags(1 ./ K, 0, ns, ns);
    res3_L          = EW_B_L * spdiags(1 ./ K, 0, ns, ns);
    res3_Bp         = EW_B_Bp * spdiags(1 ./ K, 0, ns, ns);
    res3_c_b        = EW_B_c_b * spdiags(1 ./ K, 0, ns, ns);
    res3_gam        = EW_B_gam * spdiags(1 ./ K, 0, ns, ns);

    res4_c_w        = spdiags(glob.lambda * utility_cc(c_b, 0, param, glob, options) .* c_b_c_w - ...
    ((1 - glob.lambda) * utility_cc(c_w, L, param, glob, options) - (-mu) .* B .* utility_ccc(c_w, L, param, glob, options) + (-gam_c_w) .* H_c + (-gam) .* H_c_c_w) ...
            , 0, ns, ns);
    res4_L          = spdiags(glob.lambda * utility_cc(c_b, 0, param, glob, options) .* c_b_L - ...
    ((1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - (-mu) .* B .* utility_ccl(c_w, L, param, glob, options) + (-gam_L) .* H_c + (-gam) .* H_c_L) ...
            , 0, ns, ns);
    res4_Bp         = spdiags(glob.lambda * utility_cc(c_b, 0, param, glob, options) .* c_b_Bp - ...
    ( (-gam_Bp) .* H_c + (-gam) .* H_c_Bp) ...
            , 0, ns, ns);
    res4_c_b        = spdiags(glob.lambda * utility_cc(c_b, 0, param, glob, options) .* c_b_c_b, 0, ns, ns);
    res4_gam        = spdiags(glob.lambda * utility_cc(c_b, 0, param, glob, options) .* c_b_gam - ...
    ( (-gam_gam) .* H_c), 0, ns, ns);
    
    res5_c_w        = spdiags((-gam) .* (-H_l_c_w - production_l(Z, K, L, param, glob, options) .* H_c_c_w) - ...
    ( production_l(Z, K, L, param, glob, options) .* ( (1 - glob.lambda) * utility_cc(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccc(c_w, L, param, glob, options) ) + (1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccl(c_w, L, param, glob, options) ), 0, ns, ns);  
    res5_L          = spdiags((-gam) .* (-H_l_L - production_l(Z, K, L, param, glob, options) .* H_c_L - production_ll(Z, K, L, param, glob, options) .* H_c) - ...
    ( production_l(Z, K, L, param, glob, options) .* ( (1 - glob.lambda) * utility_cl(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_ccl(c_w, L, param, glob, options) ) + ...
        production_ll(Z, K, L, param, glob, options) .* ( (1 - glob.lambda) * utility_c(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cc(c_w, L, param, glob, options) ) + (1 - glob.lambda) * utility_ll(c_w, L, param, glob, options) - ...
        (-mu) .* B .* utility_cll(c_w, L, param, glob, options) ), 0, ns, ns);
    res5_Bp         = sparse(ns, ns);
    res5_c_b        = sparse(ns, ns);
    res5_gam        = spdiags((-ones(ns, 1)) .* (-H_l - production_l(Z, K, L, param, glob, options) .* H_c), 0, ns, ns);
    
    % Put together the sparse jacobian
    jac                             = speye(3 * ns, 3 * ns);
    jac(1:ns, 1:ns)                 = res1_c_w;
    jac(1:ns, ns+1:2*ns)            = res1_L;
    jac(1:ns, 2*ns+1:3*ns)          = res1_Bp;
    jac(1:ns, 3*ns+1:4*ns)          = res1_c_b;
    jac(1:ns, 4*ns+1:5*ns)          = res1_gam;
    
    jac(ns+1:2*ns, 1:ns)            = res2_c_w;
    jac(ns+1:2*ns, ns+1:2*ns)       = res2_L;
    jac(ns+1:2*ns, 2*ns+1:3*ns)     = res2_Bp;
    jac(ns+1:2*ns, 3*ns+1:4*ns)     = res2_c_b;
    jac(ns+1:2*ns, 4*ns+1:5*ns)     = res2_gam;
    
    jac(2*ns+1:3*ns, 1:ns)          = res3_c_w;
    jac(2*ns+1:3*ns, ns+1:2*ns)     = res3_L;
    jac(2*ns+1:3*ns, 2*ns+1:3*ns)   = res3_Bp;
    jac(2*ns+1:3*ns, 3*ns+1:4*ns)   = res3_c_b;
    jac(2*ns+1:3*ns, 4*ns+1:5*ns)   = res3_gam;
    
    jac(3*ns+1:4*ns, 1:ns)          = res4_c_w;
    jac(3*ns+1:4*ns, ns+1:2*ns)     = res4_L;
    jac(3*ns+1:4*ns, 2*ns+1:3*ns)   = res4_Bp;
    jac(3*ns+1:4*ns, 3*ns+1:4*ns)   = res4_c_b;
    jac(3*ns+1:4*ns, 4*ns+1:5*ns)   = res4_gam;
    
    jac(4*ns+1:5*ns, 1:ns)          = res5_c_w;
    jac(4*ns+1:5*ns, ns+1:2*ns)     = res5_L;
    jac(4*ns+1:5*ns, 2*ns+1:3*ns)   = res5_Bp;
    jac(4*ns+1:5*ns, 3*ns+1:4*ns)   = res5_c_b;
    jac(4*ns+1:5*ns, 4*ns+1:5*ns)   = res5_gam;
    
    if strcmp(options.FB, 'Y')
        jac(1:ns, :)                = sparse(ns, 5 * ns);
        jac(2*ns+1:3*ns, :)         = sparse(ns, 5 * ns);
        jac(:, 2*ns+1:3*ns)         = sparse(5 * ns, ns);
        jac(:, 4*ns+1:5*ns)         = sparse(5 * ns, ns);
    end
end

