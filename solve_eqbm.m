function eq  = solve_eqbm(param, glob, options)

%% Globals 
s               = glob.s;  
kgrid           = glob.kgrid;
bgrid           = glob.bgrid;
ns              = size(s, 1);
Phi             = glob.Phisp;

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
for citer = 1:options.Nbackw
    % Compute current period equilibrium given the future policies
    c               = fsolve(@(x) eval_resid_backward(cold, x, param, glob, options), cold, options.optim);
    % Compute distances and update
    dc              = norm(c - cold) / norm(cold);
    cold            = c;
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
% Try to improve solving the whole system directly
cold        = fsolve(@(x) eval_resid(x, param, glob, options), cold, options.optim);

if strcmp(options.save_eqbm, 'Y')
    save(['cresult' options.GHH '.mat'], 'cold');
end

[~, ~, c_w, c_b, L, Y, I, q, Kp, Bp, r]     = eval_resid(cold, param, glob, options);     % Get eqbm objects

% if strcmp(options.print, 'Y')
%     fprintf('~~~~~ Newton iterations ~~~~~\n');
% end
% eq.flag.cconv   = false;
% for citer = 1:options.Nnewt
%     % Compute values
%     [res, jac]      = eval_resid(cold, param, glob, options);               % Evaluate residuals, get the jacobian
%     % Update c 
%     c               = cold - (jac \ res);                                       % New c
%     % Compute distances and update
%     dc              = norm(c - cold) / norm(cold);
%     cold            = c;
%     if strcmp(options.print, 'Y')
%         fprintf('%i\tdc = %1.2e\tTime: %3.2f\n', citer, dc, toc(totaltic));
%         fprintf('norm: %3.2e\n', norm(res));
%     end
%     % Check convergence
%     if (dc < options.tolc)
%         eq.flag.cconv = true;
%     end
%     if eq.flag.cconv
%         break
%     end
% end
% return;

%% Solve again on a finder grid for k
glob.Phi_Z      = glob.Phi_Zf; 

%% Compute stationary distribution
[~, ~, c_w, c_b, L, Y, I, q, Kp, Bp, r]   = eval_resid(cold, param, glob, options);
Kp              = min(Kp, glob.kmax);                                 % Restrict k's to not exceed the grid
Bp              = min(Bp, Kp * glob.bmax);                                 % Restrict k's to not exceed the grid
fspaceergk      = fundef({'spli', glob.kgrid, 0, 1},...                      % Linear approximation
                        {'spli', glob.bgrid, 0, 1});                     
QK              = funbas(fspaceergk, [Kp, Bp ./ Kp]);                                % Get basis
QZ              = glob.QZ;
Q               = dprod(QZ, QK);                                         % Full transition matrix

[vv,dd]         = eigs(Q');
dd              = diag(dd);
Lv              = vv(:,dd==max(dd));
Ldist               = Lv/sum(Lv);
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

% % Plot stationary distribution
% if strcmp(options.plotSD,'Y');
%     H = figure(888);
%     set(H,'Pos',[1          41        1920         964]);
%     Jk  = numel(glob.kgridf);
%     Jz  = numel(glob.zgridf);
%     Lz  = kron(eye(Jz),ones(1,Jk))*L;                                   % Marginals
%     Lk  = kron(ones(1,Jz),eye(Jk))*L;
%     % Marginal K
%     subplot(2,2,1);
%     plot(glob.kgridf,Lk,'o-');title('Stationary Capital Dist - Lk');grid on;
%     % Marginal Z
%     subplot(2,2,2);
%     plot(exp(glob.zgridf),Lz,'o-');title('Stationary Prod Dist - Lz');grid on;
%     eq.Lk   = Lk; 
%     eq.Lz   = Lz;
%     % Joint (K,Z) - Surface plot
%     subplot(2,2,4)
%     Lmat    = reshape(L,Jk,Jz);
%     Zmat    = repmat(glob.zgridf',Jk,1); 
%     Kmat    = repmat(glob.kgridf,1,Jz);
%     if Lk(end)<0.001;
%         Kub     = glob.kgridf(find(cumsum(Lk)>0.98,1,'first'));
%     else
%         Kub     = max(glob.kgridf);
%     end
%     Zub     = glob.zgridf(find(cumsum(Lz)>0.98,1,'first')); 
%     mesh(Zmat,Kmat,Lmat,'LineWidth',2);
%     xlabel('Productivity - z');ylabel('Capital - k');title('Joint Distribution');
%     xlim([min(glob.zgridf),Zub]);
%     ylim([min(glob.kgridf),Kub]); 
%     zlim([0,max(max(Lmat))]);
%     drawnow;
% end
% 
% %% Compute aggregates and implied price p
% Ya      = L'*v.Y;
% Ia      = L'*v.I;
% Ka      = L'*sf(:,1);
% ACa     = L'*v.AC;
% Ca      = Ya - Ia - ACa;     
% p       = 1/Ca;
% 
%% Pack-up output

eq.c_w      = c_w;
eq.c_b      = c_b;
eq.L        = L;
eq.Y        = Y;
eq.I        = I;
eq.q        = q;
eq.Kp       = Kp;
eq.Bp       = Bp;
eq.r        = r;
eq.c        = cold;
eq.dist     = Ldist;

%% Plot
c_w_arr = reshape(c_w, glob.Nk, glob.Nb, glob.Nz);
c_b_arr = reshape(c_b, glob.Nk, glob.Nb, glob.Nz);
r_arr = reshape(r, glob.Nk, glob.Nb, glob.Nz);
Bp_arr = reshape(Bp, glob.Nk, glob.Nb, glob.Nz);
Kp_arr = reshape(Kp, glob.Nk, glob.Nb, glob.Nz);
q_arr = reshape(q, glob.Nk, glob.Nb, glob.Nz);
L_arr = reshape(L, glob.Nk, glob.Nb, glob.Nz);
I_arr = reshape(I, glob.Nk, glob.Nb, glob.Nz);
w = production_l(s(:, 3), s(:, 1), L, param, glob, options);
w_arr = reshape(w, glob.Nk, glob.Nb, glob.Nz);


figure;
subplot(1, 2, 1);
plot(glob.bgrid, q_arr(8, :, 1), glob.bgrid, q_arr(9, :, 1), glob.bgrid, q_arr(10, :, 1))
xlabel('D / K')
ylabel('q')
legend(['K = ' num2str(glob.kgrid(8))], ['K = '  num2str(glob.kgrid(9))], ['K = '  num2str(glob.kgrid(10))]);
title('Capital price, Z = 0.90')

subplot(1, 2, 2);
plot(glob.bgrid, I_arr(8, :, 1), glob.bgrid, I_arr(9, :, 1), glob.bgrid, I_arr(10, :, 1))
xlabel('D / K')
ylabel('I')
legend(['K = ' num2str(glob.kgrid(8))], ['K = '  num2str(glob.kgrid(9))], ['K = '  num2str(glob.kgrid(10))]);
title('Investment, Z = 0.90')


dist_arr = reshape(Ldist, glob.Nk, glob.Nb, glob.Nz);
figure;
surf(glob.kgrid, glob.bgrid, dist_arr(:, :, 1)');
xlabel('K');
ylabel('D / K')
zlabel('Density')
title('Ergodic distribution, Z = 0.90')
end

