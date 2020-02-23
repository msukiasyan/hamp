function sim = simulate_model(s0, Nsim, T, eq, param, glob, options)

%% UNPACK
zgrid       = glob.zgrid;
c           = eq.c;
fspace      = glob.fspace;
Nz          = glob.Nz;
zub         = max(zgrid);
zlb         = min(zgrid);
rho         = glob.rhoz;
sigma       = glob.sige; 

%% SIMULATION

% Initial state
st          = s0;

% Storage
Kt          = zeros(Nsim, T);
Bt          = zeros(Nsim, T);
Zt          = zeros(Nsim, T);
ratt        = zeros(Nsim, T);
c_wt        = zeros(Nsim, T);
c_bt        = zeros(Nsim, T);
Yt          = zeros(Nsim, T);
Lt          = zeros(Nsim, T);
It          = zeros(Nsim, T);
qt          = zeros(Nsim, T);
rt          = zeros(Nsim, T);
levt        = zeros(Nsim, T);

% Initial state
Kt(:, 1)    = st(1);
Bt(:, 1)    = st(2);
Zt(:, 1)    = st(3);
ratt(:, 1)  = st(2) / st(1);

% Policy interpolant
cKp         = funfitxy(fspace, glob.s, eq.Kp); 
cBp         = funfitxy(fspace, glob.s, eq.Bp); 
cc_w        = funfitxy(fspace, glob.s, eq.c_w); 
cc_b        = funfitxy(fspace, glob.s, eq.c_b); 
cY          = funfitxy(fspace, glob.s, eq.Y); 
cL          = funfitxy(fspace, glob.s, eq.L); 
cI          = funfitxy(fspace, glob.s, eq.I); 
cq          = funfitxy(fspace, glob.s, eq.q); 
cr          = funfitxy(fspace, glob.s, eq.r); 

% Draw productivity shocks
rng(999);
PS          = randn(Nsim, T); 

% Simulation
for t = 1:T
    % Current state
    K       = Kt(:, t);
    B       = Bt(:, t);
    Z       = Zt(:, t);
    rat     = ratt(:, t);
    % Compute policies given states
    X       = funeval([cKp, cBp, cc_w, cc_b, cY, cL, cI, cq, cr], fspace, [K, rat, Z]);
    Kp      = X(:, 1); 
    Bp      = X(:, 2);
    Kp      = min(Kp, glob.kmax);
    Bp      = min(Bp, Kp .* glob.bmax);
    
    ratp    = Bp ./ Kp;
    c_w     = X(:, 3);
    c_b     = X(:, 4);
    Y       = X(:, 5);
    L       = X(:, 6);
    I       = X(:, 7);
    q       = X(:, 8);
    r       = X(:, 9);
    % Evolve Z stochastically
    Zp      = exp(rho * log(Z) + sigma * PS(:, t));
    Zp      = max(min(Zp, zub), zlb);     % Restrict to the grid
    
    %______________________________________________________________________
    % Record next period states
    if t < T
        Kt(:, t + 1)   = Kp;
        Bt(:, t + 1)   = Bp;
        Zt(:, t + 1)   = Zp;
        ratt(:, t + 1) = ratp;
    end
    %______________________________________________________________________
    % Record period t objects, these depend on...
    c_wt(:, t)     = c_w;
    c_bt(:, t)     = c_b;
    Yt(:, t)       = Y;
    Lt(:, t)       = L;
    It(:, t)       = I;
    qt(:, t)       = q;
    rt(:, t)       = r;
    levt(:, t)     = q .* K ./ (production_k(Z, K, L, param, glob, options) .* K + (1 - param.delta) * q .* K - B);
end

%% Pack
sim.Kt          = Kt;
sim.Bt          = Bt;
sim.Zt          = Zt;
sim.c_wt        = c_wt;
sim.c_bt        = c_bt;
sim.Yt          = Yt;
sim.Lt          = Lt;
sim.It          = It;
sim.qt          = qt;
sim.rt          = rt;
sim.levt        = levt;

%% Figures
% if strcmp(options.simplot,'Y')
% 
%     varlist     = {'Kt','Zt','Nt','It'};
%     titlelist   = { 'A. Capital','B. Productivity','C. Labor','D. Investment'}; 
% 
% %     h4 = figure(round(rand*1000));
% %     tvec = 0:1:T-1;
% %     for x = 1:4;
% %         subplot(2,2,x);
% %         eval(sprintf('X=%s;',char(varlist(x))));
% %         if x<=2
% %             plot(tvec,X,'-','LineWidth',2);grid on;hold on; 
% %         else
% %             plot(tvec(1:end-1),X,'-','LineWidth',2);grid on;hold on; 
% %         end
% %         xlabel('Periods','fontsize',16);
% %         title(titlelist(x),'fontsize',16);
% %         set(gca,'fontsize',16);
% %         xlim([0,T-1]);
% %     end
%     
%     varlist     = {'Kat','Cat','Nat','Iat'};
%     titlelist   = { 'A. Capital','B. Consumption','C. Labor','D. Investment'}; 
%     
%     h5 = figure(round(rand*1000));
%     tvec = 0:1:T-1;
%     for x = 1:4;
%         subplot(2,2,x);
%         eval(sprintf('X=%s;',char(varlist(x))));
%         if x<=1
%             plot(tvec,X,'-','LineWidth',2);grid on;hold on; 
%         else
%             plot(tvec(1:end-1),X(2:end),'-','LineWidth',2);grid on;hold on;  
%         end
%         xlabel('Periods','fontsize',16);
%         title(titlelist(x),'fontsize',16);
%         set(gca,'fontsize',16);
%         xlim([0,T-1]);
%     end
%     
% end
% 
