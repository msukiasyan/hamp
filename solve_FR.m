function sol  = solve_FR(eq, glob_eq, glob_old, param, glob, options)

%% Globals 
s               = glob.s;  
kgrid           = glob.kgrid;
bgrid           = glob.bgrid;
ns              = size(s, 1);
Phi             = glob.Phisp;

%% Initialise guesses
switch options.guess_FR
    case 'saved'
        load(['FR_cresult' options.GHH '.mat'], 'cold');
        if any(glob.n_FR ~= glob.n_FR_old)
            Phi_Kp          = splibas(glob_old.kgrid0, 0, glob_old.spliorder_FR(1), glob.s(:, 1));           % Basis for all k'
            Phi_Bp          = splibas(glob_old.bgrid0, 0, glob_old.spliorder_FR(2), glob.s(:, 2));           % Basis for all b'
            Phi_mup         = splibas(glob_old.mugrid0, 0, glob_old.spliorder_FR(3), glob.s(:, 3));           % Basis for all mu'
            Phi_Zp          = splibas(glob_old.zgrid0, 0, glob_old.spliorder_FR(4), glob.s(:, 4));           % Basis for all mu'
            Phi_KBmuZp      = dprod(Phi_Zp, dprod(Phi_mup, dprod(Phi_Bp, Phi_Kp)));                 % Basis for all (k', b', z')
            nso             = size(glob_old.s, 1);
            cold(1:ns)              = Phi \ (Phi_KBmuZp * cold(1:nso));
            cold(ns+1:2*ns)         = Phi \ (Phi_KBmuZp * cold(nso+1:2*nso));
            cold(2*ns+1:3*ns)       = Phi \ (Phi_KBmuZp * cold(2*nso+1:3*nso));
            cold(3*ns+1:4*ns)       = Phi \ (Phi_KBmuZp * cold(3*nso+1:4*nso));
        end
    otherwise
        c_w0            = reshape(permute(repmat(reshape(eq.c_w, glob_eq.Nk, glob_eq.Nb, glob_eq.Nz), 1, 1, 1, glob.Nmu), [1 2 4 3]), glob_eq.Nk * glob_eq.Nb * glob.Nmu * glob_eq.Nz, 1);
        c_b0            = reshape(permute(repmat(reshape(eq.c_b, glob_eq.Nk, glob_eq.Nb, glob_eq.Nz), 1, 1, 1, glob.Nmu), [1 2 4 3]), glob_eq.Nk * glob_eq.Nb * glob.Nmu * glob_eq.Nz, 1);
        L0              = reshape(permute(repmat(reshape(eq.L, glob_eq.Nk, glob_eq.Nb, glob_eq.Nz), 1, 1, 1, glob.Nmu), [1 2 4 3]), glob_eq.Nk * glob_eq.Nb * glob.Nmu * glob_eq.Nz, 1);
        Bp0             = reshape(permute(repmat(reshape(eq.Bp, glob_eq.Nk, glob_eq.Nb, glob_eq.Nz), 1, 1, 1, glob.Nmu), [1 2 4 3]), glob_eq.Nk * glob_eq.Nb * glob.Nmu * glob_eq.Nz, 1);
        s0              = gridmake({glob_eq.kgrid, glob_eq.bgrid, glob.mugrid, glob_eq.zgrid});

        meancons        = (c_w0 + c_b0) / 2;
        intp            = 0.00;
        c_w0            = c_w0 * intp + meancons * (1.0 - intp);
        c_b0            = c_b0 * intp + meancons * (1.0 - intp);

        c1old           = funfitxy(glob.fspace, s0, c_w0);
        c2old           = funfitxy(glob.fspace, s0, L0);
        c3old           = funfitxy(glob.fspace, s0, Bp0);
        c4old           = funfitxy(glob.fspace, s0, c_b0);

        cold            = [c1old; c2old; c3old; c4old];
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
coldold = cold;
for citer = 1:options.Nbackw
    % Compute current period equilibrium given the future policies
    cnum                    = zeros(4*ns, 1);
    c                       = zeros(4*ns, 1);
    cnum(1:ns)              = Phi * cold(1:ns);
    cnum(ns+1:2*ns)         = Phi * cold(ns+1:2*ns);
    cnum(2*ns+1:3*ns)       = Phi * cold(2*ns+1:3*ns);
    cnum(3*ns+1:4*ns)       = Phi * cold(3*ns+1:4*ns);
    cnumnew               = fsolve(@(x) eval_resid_backward_FR(cold, x, param, glob, options), cnum, optim_fast);
    c(1:ns)              = Phi \ cnumnew(1:ns);
    c(ns+1:2*ns)         = Phi \ cnumnew(ns+1:2*ns);
    c(2*ns+1:3*ns)       = Phi \ cnumnew(2*ns+1:3*ns);
    c(3*ns+1:4*ns)       = Phi \ cnumnew(3*ns+1:4*ns);
    % Compute distances and update
    dc              = norm(c - cold) / norm(cold);
    cold            = c;
    if strcmp(options.print, 'Y')
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n', citer, dc, toc(totaltic));
        fprintf('norm: %3.2e\n', norm(eval_resid_FR(cold, param, glob, options)));
    end
    % Check convergence
    if (dc < options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end
% Try to improve solving the whole system directly
cold        = fsolve(@(x) eval_resid_FR(x, param, glob, options), cold, options.optim);

if strcmp(options.save_eqbm_FR, 'Y')
    save(['FR_cresult' options.GHH '.mat'], 'cold');
    glob.n_FR_old           = glob.n_FR;
end

[~, ~, ~, ~, c_w, c_b, L, Y, I, q, Kp, Bp]     = eval_resid_FR(cold, param, glob, options);     % Get eqbm objects


%% Solve again on a finer grid for k
% glob.Phi_Z      = glob.Phi_Zf; 

%% Compute stationary distribution
[~, ~, ~, ~, c_w, c_b, L, Y, I, q, Kp, Bp, gam]   = eval_resid_FR(cold, param, glob, options);
Kp              = min(Kp, glob.kmax);                                 % Restrict k's to not exceed the grid
Bp              = min(Bp, Kp * glob.bmax);                                 % Restrict k's to not exceed the grid
mup             = min(gam, glob.mumax_FR);                                 % 
fspaceergk      = fundef({'spli', glob.kgrid, 0, 1},...                      % Linear approximation
                        {'spli', glob.bgrid, 0, 1}, ...
                        {'spli', glob.mugrid, 0, 1});                     
QK              = funbas(fspaceergk, [Kp, Bp ./ Kp, mup]);                                % Get basis
QZ              = glob.QZ;
Q               = dprod(QZ, QK);                                         % Full transition matrix

[vv,dd]         = eigs(Q');
dd              = diag(dd);
Lv              = vv(:,dd==max(dd));
Ldist           = Lv/sum(Lv);
% Ldist           = ones(size(Q, 1), 1);                                    % Start from uniform dist
% Ldist           = Ldist / sum(Ldist);
% 
% for itL = (1:options.itermaxL)                                         % Iterate to find stat dist
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
 
%% Pack-up output

sol.c_w      = c_w;
sol.c_b      = c_b;
sol.L        = L;
sol.Y        = Y;
sol.I        = I;
sol.q        = q;
sol.Kp       = Kp;
sol.Bp       = Bp;
sol.c        = cold;
sol.dist     = Ldist;

%% Plot
c_w_arr     = reshape(c_w, glob.Nk, glob.Nb, glob.Nz);
c_b_arr     = reshape(c_b, glob.Nk, glob.Nb, glob.Nz);
r_arr       = reshape(r, glob.Nk, glob.Nb, glob.Nz);
Bp_arr      = reshape(Bp, glob.Nk, glob.Nb, glob.Nz);
Kp_arr      = reshape(Kp, glob.Nk, glob.Nb, glob.Nz);
q_arr       = reshape(q, glob.Nk, glob.Nb, glob.Nz);
L_arr       = reshape(L, glob.Nk, glob.Nb, glob.Nz);
I_arr       = reshape(I, glob.Nk, glob.Nb, glob.Nz);
w           = production_l(s(:, 3), s(:, 1), L, param, glob, options);
w_arr       = reshape(w, glob.Nk, glob.Nb, glob.Nz);


figure;
subplot(1, 2, 1);
plot(glob.bgrid, q_arr(8, :, 1), glob.bgrid, q_arr(9, :, 1), glob.bgrid, q_arr(10, :, 1), 'LineWidth', 2)
xlabel('D / K')
ylabel('q')
legend(['K = ' num2str(glob.kgrid(8))], ['K = '  num2str(glob.kgrid(9))], ['K = '  num2str(glob.kgrid(10))]);
title('Capital price, Z = 0.90')

subplot(1, 2, 2);
plot(glob.bgrid, I_arr(8, :, 1), glob.bgrid, I_arr(9, :, 1), glob.bgrid, I_arr(10, :, 1), 'LineWidth', 2)
xlabel('D / K')
ylabel('I')
legend(['K = ' num2str(glob.kgrid(8))], ['K = '  num2str(glob.kgrid(9))], ['K = '  num2str(glob.kgrid(10))]);
title('Investment, Z = 0.90')

figure;
subplot(1, 2, 1);
plot(glob.bgrid, L_arr(8, :, 1), glob.bgrid, L_arr(9, :, 1), glob.bgrid, L_arr(10, :, 1), 'LineWidth', 2)
xlabel('D / K')
ylabel('L')
legend(['K = ' num2str(glob.kgrid(8))], ['K = '  num2str(glob.kgrid(9))], ['K = '  num2str(glob.kgrid(10))]);
title('Labor, Z = 0.90')

subplot(1, 2, 2);
plot(glob.bgrid, r_arr(8, :, 1), glob.bgrid, r_arr(9, :, 1), glob.bgrid, r_arr(10, :, 1), 'LineWidth', 2)
xlabel('D / K')
ylabel('r')
legend(['K = ' num2str(glob.kgrid(8))], ['K = '  num2str(glob.kgrid(9))], ['K = '  num2str(glob.kgrid(10))]);
title('Risk-free rate, Z = 0.90')


dist_arr = reshape(Ldist, glob.Nk, glob.Nb, glob.Nz);
figure;
surf(glob.kgrid, glob.bgrid, dist_arr(:, :, 1)');
xlabel('K');
ylabel('D / K')
zlabel('Density')
title('Ergodic distribution, Z = 0.90')
end

