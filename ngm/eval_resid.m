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
Z               = glob.s(:, 2);
% Unpack
c1              = c(1:ns);
c               = Phi * c1;

c               = max(c, 0);

%% Solve equations 
L               = solve_L_from_MRS(Z, K, c, param, glob, options);
Y               = production(Z, K, L, param, glob, options);
I               = Y - c;
q               = 1 ./ cap_prod_prime(I ./ K, param, glob, options);
Kp              = K .* (cap_prod(I ./ K, param, glob, options) + 1 - param.delta);
mpk             = production_k(Z, K, L, param, glob, options);
mu              = utility_c(c, L, param, glob, options);
mu_prod         = mu .* (mpk + (1 - param.delta) * q);

% Enforce bounds
% Kp              = min(max(Kp, glob. kmin), glob.kmax);
% Bp              = min(max(Bp, glob. bmin), glob.bmax);

%% Create basis matrices for next states
Phi_Kp          = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp);           % Basis for all k'
Phi_KZp         = dprod(glob.Phi_Z, Phi_Kp);                 % Basis for all (k', b', z')

%% Compute next state objects
PhiEmu_prod     = basiscast * mu_prod;                                     % Approximate on the same basis, take expectation and get coefficients again                    
Emu_prodp       = Phi_KZp * PhiEmu_prod;                               % Evaluate at (k', b', z)


%% Compute residuals
res                 = zeros(1 * ns, 1);
res(1:ns)           = q .* c - utility_c_inv(glob.beta * Emu_prodp, L, param, glob, options);                     % Euler equation for workers

%% Compute the jacobian if requested
jac             = [];
if (nargout == 2)
    % Compute derivatives with respect to all three inputs
    L_c             = solve_L_from_MRS_c(Z, K, c, param, glob, options);
    
    Y_c             = production_l(Z, K, L, param, glob, options) .* L_c; % Chain rule
    
    I_c             = Y_c - ones(ns, 1);
    
    q_c             = -q .^ 2 .* cap_prod_prime_prime(I ./ K, param, glob, options) ./ K .* I_c;
    
    Kp_c            = (1 ./ q) .* I_c;
    
    mpk_c           = production_kl(Z, K, L, param, glob, options) .* L_c;
    
    mu_c            = utility_cc(c, L, param, glob, options) .* ones(ns, 1) + utility_cl(c, L, param, glob, options) .* L_c;
    
    mu_prod_c       = mu_c .* (mpk + (1 - param.delta) * q) + mu .* (mpk_c + (1 - param.delta) * q_c);
    
    % Create derivative matrices
    Phi_Kp_der      = splibas(glob.kgrid0, 0, glob.spliorder(1), Kp, 1);
    Phi_KZp_Kder    = dprod(glob.Phi_Z, Phi_Kp_der);       
    
    Kder_prod       = spdiags(Phi_KZp_Kder * PhiEmu_prod, 0, ns, ns);
    
    diagKp_c        = spdiags(Kp_c, 0, ns, ns);
    
    next_basis      = Phi_KZp * basiscast;
   

    Emu_prodp_c     = Kder_prod * diagKp_c + ... 
                        next_basis * spdiags(mu_prod_c, 0, ns, ns);

    % Calculate the blocks of the jacobian

    mu_prime1       = spdiags(utility_c_inv_prime(glob.beta * Emu_prodp, L, param, glob, options), 0, ns, ns);

    res1_c          = spdiags(q, 0, ns, ns) + spdiags(q_c .* c, 0, ns, ns) - mu_prime1 * glob.beta * Emu_prodp_c;
    
    % Put together the jacobian
    jac                             = zeros(ns, ns);
    jac(1:ns, 1:ns)                 = res1_c * Phi;
end

