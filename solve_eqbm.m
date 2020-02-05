function eq  = solve_eqbm(param, glob, options)

%% Globals 
s               = glob.s;  
kgrid           = glob.kgrid;
bgrid           = glob.kgrid;
ns              = size(s, 1);
Phi             = glob.Phisp;

%% Initialise guesses (if val.cresult has an old guess in it, use that)
c_w0            = ones(ns, 1) * 5; % ones(ns, 1) / 2;                                             % Initial guess for c_w
c_b0            = ones(ns, 1) * 5; % ones(ns, 1) / 2;                                             % Initial guess for c_b
r0              = ones(ns, 1);                                              % Initial guess for r
c1old           = Phi \ c_w0;                                               % Coefficients for c_w
c2old           = Phi \ c_b0;                                               % Coefficients for c_b
c3old           = Phi \ r0;                                               % Coefficients for c_b
if isfield(options, 'cresult') && ~isempty(options.cresult)
    c1old       = options.cresult(1:ns);
    c2old       = options.cresult(ns + 1:2 * ns);
    c3old       = options.cresult(2 * ns + 1:end);
end
cold            = [c1old; c2old; c3old];
totaltic        = tic;

%% Newton iterations
% c               = fsolve(@(x) eval_resid(x, param, glob, options), cold, optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true));
% return

if strcmp(options.print, 'Y')
    fprintf('~~~~~ Newton iterations ~~~~~\n');
end
eq.flag.cconv   = false;
for citer = 1:options.Nnewt
    % Compute values
    [res, jac]      = eval_resid(cold, param, glob, options);               % Evaluate residuals, get the jacobian
    % Update c 
    c               = cold - 0.01 * (jac \ res);                                       % New c
    % Compute distances and update
    dc              = norm(c - cold) / norm(cold);
    cold            = c;
    if strcmp(options.print, 'Y')
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n', citer, dc, toc(totaltic));
        fprintf('norm: %3.2f\n', norm(res));
    end
    % Check convergence
    if (dc < options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end
return;
% %% Solve again on a finder grid for k
% glob.Phi_Z      = glob.Phi_Zf; 
% v               = solve_valfunc(c,sf,P,param,glob,options,1);           % Solve using the obtained approximation but on many points
% 
% %% Compute stationary distribution
% Kp              = min(v.Kp,max(kgrid));                                 % Restrict k's to not exceed the grid
% fspaceergk      = fundef({'spli',glob.kgridf,0,1});                     % Linear approximation
% QK              = funbas(fspaceergk,Kp);                                % Get basis
% QZ              = glob.QZ;
% Q               = dprod(QZ,QK);                                         % Full transition matrix
% 
% % [vv,dd]         = eigs(Q');
% % dd              = diag(dd);
% % Lv              = vv(:,dd==max(dd));
% % L               = Lv/sum(Lv);
% L               = ones(size(Q,1),1);                                    % Start from uniform dist
% L               = L/sum(L);
% 
% for itL = (1:options.itermaxL);                                         % Iterate to find stat dist
%     Lnew    = Q'*L;  
%     dL      = norm(Lnew-L)/norm(L);  
%     if (dL<options.tolL),break,end;
%     if mod(itL,100)==0 
%         if strcmp(options.print,'Y')
%             fprintf('dL:\t%1.3e\n',dL);
%         end
%     end
%     L       = Lnew;
% end
% 
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
% %% Pack-up output
% eq.v    = v;
% eq.c    = c;
% eq.p    = p;
% eq.L    = L;
% eq.Ya   = Ya;
% eq.Ca   = Ca;
% eq.Ia   = Ia;
% eq.ACa  = ACa;
% eq.Ka   = Ka;
% eq.Q    = Q;
% 
% end
% 
