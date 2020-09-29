function eq  = solve_eqbm(param, glob, options)

%% Globals 
s               = glob.s;  
kgrid           = glob.kgrid;
bgrid           = glob.bgrid;
ns              = size(s, 1);
Phi             = glob.Phi;
Phil            = glob.Phil;
Phiu            = glob.Phiu;
Phif            = glob.Phif;

%% Initialise guesses
switch options.guess
    case 'file'
        load('init_guess.mat', 'fres');

        fres            = fres(fres(:, 2) > 0, :);
        Nf              = size(fres, 1);

        k0              = kron(ones(glob.Nz, 1), fres(:, 1));
        b0              = kron(ones(glob.Nz, 1), fres(:, 2));
        c_w0            = kron(ones(glob.Nz, 1), fres(:, 3));
        c_b0            = kron(ones(glob.Nz, 1), fres(:, 4));
        r0              = kron(ones(glob.Nz, 1), fres(:, 5));
        z0              = kron(glob.zgrid, ones(Nf, 1));
        % c_b0(c_b0 < 1e-5)   = 1e-5;

        c1old           = funfitxy(glob.fspace, [k0, b0 ./ k0, z0], c_w0);
        c2old           = funfitxy(glob.fspace, [k0, b0 ./ k0, z0], c_b0);
        c3old           = funfitxy(glob.fspace, [k0, b0 ./ k0, z0], r0);
        cold            = [c1old; c2old; c3old];
    case 'saved'
        load(['cresult' options.GHH '.mat'], 'cold');
    otherwise
        if options.GHH == 'N'
            c_w0            = (s(:, 2) .* s(:, 1)) * 0.1 + s(:, 1) * 0.02;                                % Initial guess for c_w
        else
            l0              = solve_L_from_MRS(s(:, 3), s(:, 1), 0, param, glob, options);
            c_w0            = (s(:, 2) .* s(:, 1)) * 0.1 + production(s(:, 3), s(:, 1), l0, param, glob, options);                                % Initial guess for c_w
        end
        c_b0            = (s(:, 1) - s(:, 2) .* s(:, 1)) * 0.02 + 0.01;                                % Initial guess for c_b
        r0              = ones(ns, 1);                                      % Initial guess for r
        c1old           = Phi \ c_w0;                                       % Coefficients for c_w
        c2old           = Phi \ c_b0;                                       % Coefficients for c_b
        c3old           = Phi \ r0;                                         % Coefficients for c_b
        cold            = [c1old; c2old; c3old];
end

totaltic        = tic;

%% Backward iterations
if strcmp(options.print, 'Y')
    fprintf('~~~~~ Backward iterations ~~~~~\n');
end
eq.flag.cconv   = false;
cnum        = [Phi * cold(1:ns); Phi * cold(ns+1:2*ns); Phi * cold(2*ns+1:3*ns)];
for citer = 1:options.Nbackw
    % Compute current period equilibrium given the future policies
    cnumnew         = fsolve(@(x) eval_resid_backward(cnum, x, param, glob, options), cnum, options.optim);
    if options.lev_cap == 'Y'
        [~, ~, ~, ~, ~, ~, ~, q, Kp, Bp, ~]     = eval_resid_backward(cnum, cnumnew, param, glob, options);
        mask        = Bp - param.chi_lev * q .* Kp >= 0;
        if sum(mask) > 0
            c_wm        = cnumnew(1:ns) * 1.1;
            c_bm        = cnumnew(ns+1:2*ns) * 1.2;
            rm          = cnumnew(2*ns+1:3*ns) * 0.95;
            cmask       = fsolve(@(x) eval_resid_backward(cnum, x, param, glob, options, glob.s(mask, :), 'Y', mask), [c_wm(mask);c_bm(mask);rm(mask)], options.optim);
            cnumnew([mask;mask;mask;])      = cmask;
        end
    end
    % Compute distances and update
    dc              = norm(cnumnew - cnum) / norm(cnum);
    cnum            = cnumnew;
    if strcmp(options.print, 'Y')
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n', citer, dc, toc(totaltic));
        fprintf('norm: %3.2e\n', norm(eval_resid(cold, param, glob, options)));
    end
    % Check convergence
    if (dc < options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end
cold        = [Phi \ cnum(1:ns); Phi \ cnum(ns+1:2*ns); Phi \ cnum(2*ns+1:3*ns)];
% Try to improve solving the whole system directly
cold        = fsolve(@(x) eval_resid(x, param, glob, options), cold, options.optim);

if strcmp(options.save_eqbm, 'Y')
    save(['cresult' options.GHH '.mat'], 'cold');
end

[~, ~, c_w, c_b, L, Y, I, q, Kp, Bp, r]     = eval_resid(cold, param, glob, options);     % Get eqbm objects

%% Compute extra variables
net_worth       = q .* Kp - Bp ./ r + c_b;
lev             = q .* Kp ./ net_worth;


%% Solve again on a finer grid for k
% glob.Phi_Z      = glob.Phi_Zf; 

%% Compute stationary distribution
[~, ~, c_w, c_b, L, Y, I, q, Kp, Bp, r]   = eval_resid(cold, param, glob, options);
Kp              = min(Kp, glob.kmax);                                 % Restrict k's to not exceed the grid
Bp              = max(min(Bp, Kp * glob.bmax), glob.bmin);

% Approximate on the finer grid
cKp             = funfitxy(glob.fspace, glob.s, Kp);
cBp             = funfitxy(glob.fspace, glob.s, Bp);
X               = funeval([cKp, cBp], glob.fspace, glob.sf);
Kpf             = X(:, 1);
Bpf             = X(:, 2);

fspaceergk      = fundef({'spli', glob.kgridf, 0, 1},...                      % Linear approximation
                        {'spli', glob.bgridf, 0, 1});
QK              = funbas(fspaceergk, [Kpf, Bpf ./ Kpf]);
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

%% Approximate solution objects on the finer grid in order to compute stats over the ergodic distribution
cc_w            = Phil \ (Phiu \ c_w);
cc_b            = Phil \ (Phiu \ c_b);
cL              = Phil \ (Phiu \ L);
cY              = Phil \ (Phiu \ Y);
cI              = Phil \ (Phiu \ I);
cq              = Phil \ (Phiu \ q);
cr              = Phil \ (Phiu \ r);

c_wf            = Phif * cc_w;
c_bf            = Phif * cc_b;
Lf              = Phif * cL;
Yf              = Phif * cY;
If              = Phif * cI;
qf              = Phif * cq;
rf              = Phif * cr;
 
%% Pack-up output

eq.c_w          = c_w;
eq.c_b          = c_b;
eq.L            = L;
eq.Y            = Y;
eq.I            = I;
eq.q            = q;
eq.Kp           = Kp;
eq.Bp           = Bp;
eq.r            = r;
eq.c            = cold;
eq.dist         = Ldist;
eq.net_worth    = net_worth;
eq.lev          = lev;
eq.Kpf          = Kpf;
eq.Bpf          = Bpf;
eq.c_wf         = c_wf;
eq.c_bf         = c_bf;
eq.Lf           = Lf;
eq.Yf           = Yf;
eq.If           = If;
eq.qf           = qf;
eq.rf           = rf;

end

