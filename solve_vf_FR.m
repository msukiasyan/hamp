function sol  = solve_vf_FR(eq, glob_eq, glob_old, param, glob, options)

%% Globals 
s               = glob.s;  
kgrid           = glob.kgrid;
bgrid           = glob.bgrid;
ns              = size(s, 1);
Phi             = glob.Phi;
Phil            = glob.Phil;
Phiu            = glob.Phiu;
Phif            = glob.Phif;
Emat            = glob.Emat;

K               = glob.s(:, 1);
B               = glob.s(:, 1) .* glob.s(:, 2);
mu              = glob.s(:, 3);
Z               = glob.s(:, 4);

%% Initialise guesses 
switch options.guess_FR
    case 'saved'
        if strcmp(options.FB, 'Y')
            load(['FR_cresult' options.GHH '_FB.mat'], 'cold', 'ecold', 'eq');
        else
            load(['FR_cresult' options.GHH '_vf.mat'], 'cold', 'ecold', 'eq');
        end
        eq0.c_w         = eq.c_w;
        eq0.L           = eq.L;
        eq0.Bp          = eq.Bp;
        eq0.c_b         = eq.c_b;
        eq0.mup         = eq.mup;
        % If the grid has changed, approximate guesses at gridpoints using the old grid
        if glob.mumax_FR ~= glob_old.mumax_FR || glob.bmax_FR ~= glob_old.bmax_FR || ...
                 glob.Ns ~= glob_old.Ns || any(any(glob.s ~= glob_old.s))
            Phi_Kp          = splibas(glob_old.kgrid0, 0, glob_old.spliorder_FR(1), glob.s(:, 1));              % Basis for all k
            Phi_Bp          = splibas(glob_old.bgrid0, 0, glob_old.spliorder_FR(2), glob.s(:, 2));              % Basis for all b
            Phi_mup         = splibas(glob_old.mugrid0, 0, glob_old.spliorder_FR(3), glob.s(:, 3));             % Basis for all mu
            Phi_Zp          = splibas(glob_old.zgrid0, 0, glob_old.spliorder_FR(4), glob.s(:, 4));              % Basis for all z
            Phi_KBmuZp      = dprod(Phi_Zp, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp)));                             % Basis for all (k, b, mu, z)
            
            cold            = Phi \ (Phi_KBmuZp * cold);
            ecold           = Phi \ (Emat * cold);
            eq0.c_w         = Phi_KBmuZp * (glob_old.Phi \ eq.c_w);
            eq0.L           = Phi_KBmuZp * (glob_old.Phi \ eq.L);
            eq0.Bp          = Phi_KBmuZp * (glob_old.Phi \ eq.Bp);
            eq0.c_b         = Phi_KBmuZp * (glob_old.Phi \ eq.c_b);
            eq0.mup         = Phi_KBmuZp * (glob_old.Phi \ eq.mup);
        end
    otherwise
        % Use the competitive eqbm values as guesses; assumes high inverse Frisch
        Phi_Kp          = splibas(glob_eq.kgrid0, 0, glob_eq.spliorder(1), glob.s(:, 1));           % Basis for all k
        Phi_Bp          = splibas(glob_eq.bgrid0, 0, glob_eq.spliorder(2), glob.s(:, 2));           % Basis for all b
        Phi_Zp          = splibas(glob_eq.zgrid0, 0, glob_eq.spliorder(3), glob.s(:, 4));           % Basis for all mu
        Phi_KBmuZp      = dprod(Phi_Zp, dprod(Phi_Bp, Phi_Kp));                                     % Basis for all (k, b, z)
        c_w0            = Phi_KBmuZp * (glob_eq.Phi \ eq.c_w);
        c_b0            = Phi_KBmuZp * (glob_eq.Phi \ eq.c_b);
        L0              = Phi_KBmuZp * (glob_eq.Phi \ eq.L);

        meancons        = (c_w0 + c_b0) / 2;
        intp            = 0.0;
        c_w0            = c_w0 * intp + meancons * (1.0 - intp);
        c_b0            = c_b0 * intp + meancons * (1.0 - intp);
        
        eq0.c_w         = c_w0;
        eq0.c_b         = c_b0;
        eq0.L           = ones(ns, 1);
        eq0.Kp          = Phi_KBmuZp * (glob_eq.Phi \ eq.Kp);
        eq0.Bp          = Phi_KBmuZp * (glob_eq.Phi \ eq.Bp) * 0.85 + (glob.s(:, 1) .* glob.s(:, 2)) * 0.15;
        eq0.mup         = glob.s(:, 3) * 0.0;
        
        c               = compute_planner_vf(eq0, param, glob, options);
        cold            = c(1:ns);
        ecold           = c(ns+1:2*ns);
        
        % For the first iteration, calculate c_b and gam exactly given the other guesses.
        % Turn these into residual equations for later iterations since
        % utility_c_inv will be undefined for negative arguments
        H_l             = (B - eq0.c_w) .* utility_cl(eq0.c_w, eq0.L, param, glob, options) -  ...
            utility_ll(eq0.c_w, eq0.L, param, glob, options) .* eq0.L - utility_l(eq0.c_w, eq0.L, param, glob, options);
        H_c             = (B - eq0.c_w) .* utility_cc(eq0.c_w, eq0.L, param, glob, options) -  ...
            utility_c(eq0.c_w, eq0.L, param, glob, options) - utility_cl(eq0.c_w, eq0.L, param, glob, options) .* eq0.L;
        gam             = -((production_l(Z, K, eq0.L, param, glob, options) .* ((1 - glob.lambda) * utility_c(eq0.c_w, eq0.L, param, glob, options) - ...
    (-mu) .* B .* utility_cc(eq0.c_w, eq0.L, param, glob, options)) + (1 - glob.lambda) * utility_l(eq0.c_w, eq0.L, param, glob, options) - ...
    (-mu) .* B .* utility_cl(eq0.c_w, eq0.L, param, glob, options)) ./ (-H_l - production_l(Z, K, eq0.L, param, glob, options) .* H_c));
        eq0.mup         = gam;
        eq0.c_b         = utility_c_inv(max(((1 - glob.lambda) * utility_c(eq0.c_w, eq0.L, param, glob, options) - ...
     (-mu) .* B .* utility_cc(eq0.c_w, eq0.L, param, glob, options) + (-gam) .* H_c) / glob.lambda, 0), 0, param, glob, options);         
end

totaltic        = tic;

%% Backward iterations
optim       = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'trust-region-dogleg', ...
    'SpecifyObjectiveGradient',true, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-30, ...
    'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIterations', 2000);
optim_fast  = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'trust-region', ... 
    'SpecifyObjectiveGradient',true, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-30, ...
    'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIterations', 2000);

if strcmp(options.print, 'Y')
    fprintf('~~~~~ Backward iterations ~~~~~\n');
end
eq.flag.cconv   = false;
coldold         = cold;
nuseq           = [(param.nu + linspace(0, 1, param.Nseq) .^ (1 / param.nucurv) * (param.nufinal - param.nu)),...
            repmat(param.nufinal, 1, options.Nbackw - param.Nseq)];
for citer = 1:options.Nbackw
    param.nu        = nuseq(citer);
    % Compute current period equilibrium given the future policies
    cnum                    = zeros(5*ns, 1);
    cnum(1:ns)              = eq0.c_w;
    cnum(ns+1:2*ns)         = eq0.L;
    cnum(2*ns+1:3*ns)       = eq0.Bp;
    cnum(3*ns+1:4*ns)       = eq0.c_b;
    cnum(4*ns+1:5*ns)       = eq0.mup;
    cnumnew                 = fsolve(@(x) eval_vf_resid_backward_FR(ecold, x, param, glob, options), cnum, optim_fast);
    [~, ~, eq1]             = eval_vf_resid_backward_FR(ecold, cnumnew, param, glob, options);
    eq0                     = eq1;
    if citer <= param.Nseq + 50
        c                       = update_planner_vf(eq1, ecold, param, glob, options);
    else                                                                    % switch to Howard improvement after enough iterations
        c                       = compute_planner_vf(eq1, param, glob, options);
        c                       = c(1:ns);
    end
    ec                      = Phi \ (Emat * c);

    % Compute distances and update
    dc              = norm(Phi * c - Phi * cold) / norm(Phi * cold) * 1e2;
    cold            = c;
    ecold           = ec;
    if strcmp(options.print, 'Y')
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n', citer, dc, toc(totaltic));
    end
    % Check convergence
    if dc < options.tolc * 1e-3 || (citer > param.Nseq && dc < options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end
c               = [cold; ecold];
% Try to improve solving the functional equation directly
c               = fsolve(@(x) eval_vf_resid(x, cnumnew, param, glob, options), c, optim_fast);
cold            = c(1:ns);
ecold           = c(ns+1:2*ns);

[~, ~, eq]      = eval_vf_resid(c, cnumnew, param, glob, options);

if strcmp(options.save_eqbm_FR, 'Y')
    if strcmp(options.FB, 'Y')
        save(['FR_cresult' options.GHH '_FB.mat'], 'cold', 'ecold', 'eq');
    else
        save(['FR_cresult' options.GHH '_vf.mat'], 'cold', 'ecold', 'eq');
    end
end

eq              = calc_implied_taxes(eq, param, glob, options);

eq.cold         = cold;
eq.ecold        = ecold;

%% Compute other stats
% Unconditional log consumption variance
c_w_arr         = reshape(eq.c_w, glob.Nk, glob.Nb, glob.Nmu, glob.Nz);
c_b_arr         = reshape(eq.c_b, glob.Nk, glob.Nb, glob.Nmu, glob.Nz);

eq.c_w_var      = std(log(c_w_arr), glob.Pssz, 4) .^ 2;
eq.c_b_var      = std(log(c_b_arr), glob.Pssz, 4) .^ 2;

%% Compute stationary distribution
Kp              = eq.Kp;
Bp              = eq.Bp;
mup             = eq.mup;
% Restrict states to not exceed the grids
Kp              = min(Kp, glob.kmax_FR);                                    
Bp              = max(min(Bp, Kp * glob.bmax_FR), glob.bmin_FR);
mup             = max(min(mup, glob.mumax_FR), glob.mumin_FR);

% Approximate on the finer grid
cKp             = Phil \ (Phiu \ Kp);
cBp             = Phil \ (Phiu \ Bp);
cmup            = Phil \ (Phiu \ mup);
Kpf             = Phif * cKp;
Bpf             = Phif * cBp;
mupf            = Phif * cmup;

fspaceergk      = fundef({'spli', glob.kgridf, 0, 1},...                        % Linear approximation
                        {'spli', glob.bgridf, 0, 1},...
                        {'spli', glob.mugridf, 0, 1});
QK              = funbas(fspaceergk, [Kpf, Bpf ./ Kpf, mupf]);
QZ              = glob.QZ;
Q               = dprod(QZ, QK);                                                % Full transition matrix

[vv,dd]         = eigs(Q');
dd              = diag(dd);
Lv              = vv(:,dd==max(dd));
Ldist           = Lv/sum(Lv);

% Ldist           = ones(size(Q, 1), 1);                                          % Start from uniform dist
% Ldist           = Ldist / sum(Ldist);
% 
% for itL = (1:options.itermaxL)                                                  % Iterate to find stat dist
%     Lnew    = Q' * Ldist;  
%     dL      = norm(Lnew - Ldist) / norm(Ldist);  
%     if (dL < options.tolL),break,end
%     if mod(itL, 100)==0 
%         if strcmp(options.print, 'Y')
%             fprintf('dL:\t%1.3e\n', dL);
%         end
%     end
%     Ldist   = Lnew;
% end

%% Approximate solution objects on the finer grid in order to compute stats over the ergodic distribution
cc_w            = Phil \ (Phiu \ eq.c_w);
cc_b            = Phil \ (Phiu \ eq.c_b);
cL              = Phil \ (Phiu \ eq.L);
cY              = Phil \ (Phiu \ eq.Y);
cI              = Phil \ (Phiu \ eq.I);
cq              = Phil \ (Phiu \ eq.q);
cr              = Phil \ (Phiu \ eq.r);
ctax_lab        = Phil \ (Phiu \ eq.implied_tax_lab);
ctax_r          = Phil \ (Phiu \ eq.implied_tax_r);
ctransfer       = Phil \ (Phiu \ eq.transfer);

c_wf            = Phif * cc_w;
c_bf            = Phif * cc_b;
Lf              = Phif * cL;
Yf              = Phif * cY;
If              = Phif * cI;
qf              = Phif * cq;
rf              = Phif * cr;
tax_labf        = Phif * ctax_lab;
tax_rf          = Phif * ctax_r;
transferf       = Phif * ctransfer;
 
%% Pack-up output
sol             = eq;
sol.Kpf         = Kpf;
sol.Bpf         = Bpf;
sol.mupf        = mupf;
sol.c_wf        = c_wf;
sol.c_bf        = c_bf;
sol.Lf          = Lf;
sol.Yf          = Yf;
sol.If          = If;
sol.qf          = qf;
sol.rf          = rf;
sol.cold        = cold;
sol.ecold       = ecold;
sol.dist        = Ldist;
sol.tax_labf    = tax_labf;
sol.tax_rf      = tax_rf;
sol.transferf   = transferf;

end

