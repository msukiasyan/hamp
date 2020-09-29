function girf = girf_planner(s0, sol, param, glob, options)
%% Unpack
Nsim        = options.Ngirf;
T           = options.Tgirf;
fspace      = glob.fspace;
Phi         = glob.Phi;
Phil        = glob.Phil;
Phiu        = glob.Phiu;
zgrid       = glob.zgrid;
P           = glob.P;

%% Simulate

% Initial state
st          = s0;

% Storage
Kt          = zeros(Nsim, T);
Bt          = zeros(Nsim, T);
mut         = zeros(Nsim, T);
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
transft     = zeros(Nsim, T);
tax_labt    = zeros(Nsim, T);
tax_rt      = zeros(Nsim, T);

% Initial state
Kt(:, 1)    = st(1);
Bt(:, 1)    = st(2);
mut(:, 1)   = st(3);
ratt(:, 1)  = st(2) / st(1);

% Policy interpolants
cKp             = Phil \ (Phiu \ sol.Kp);
cBp             = Phil \ (Phiu \ sol.Bp);
cmup            = Phil \ (Phiu \ sol.mup);
cc_w            = Phil \ (Phiu \ sol.c_w);
cc_b            = Phil \ (Phiu \ sol.c_b);
cL              = Phil \ (Phiu \ sol.L);
cY              = Phil \ (Phiu \ sol.Y);
cI              = Phil \ (Phiu \ sol.I);
cq              = Phil \ (Phiu \ sol.q);
cr              = Phil \ (Phiu \ sol.r);
ctax_lab        = Phil \ (Phiu \ sol.implied_tax_lab);
ctax_r          = Phil \ (Phiu \ sol.implied_tax_r);
ctransf         = Phil \ (Phiu \ sol.transfer);

% Draw paths from the Markov chain
rng(999);
Zi          = find(zgrid == s0(4));
simcnt      = zeros(glob.Nz, 1);
simcnt(Zi)  = Nsim;
Zi          = simulate(dtmc(P), T - 1, 'X0', simcnt)';
Zt          = zgrid(Zi);

X0          = funeval([cKp, cBp, cmup, cc_w, cc_b, cY, cL, cI, cq, cr, ctax_lab, ctax_r, ctransf], fspace, [st(1), st(2) / st(1), st(3), median(zgrid)]);

% Simulation
for t = 1:T
    % Current state
    K       = Kt(:, t);
    B       = Bt(:, t);
    mu      = mut(:, t);
    Z       = Zt(:, t);
    rat     = ratt(:, t);
    % Compute policies given states
    X       = funeval([cKp, cBp, cmup, cc_w, cc_b, cY, cL, cI, cq, cr, ctax_lab, ctax_r, ctransf], fspace, [K, rat, mu, Z]);
    
    Kp          = X(:, 1);
    Bp          = X(:, 2);
    mup         = X(:, 3);
    
    ratp        = Bp ./ Kp;
%     qual_shock  = (glob.dis_qual * (Zi == 1) + 1 * (Zi ~= 1));
%     K           = K .* qual_shock;
    
    c_w         = X(:, 4);
    c_b         = X(:, 5);
    Y           = X(:, 6);
    L           = X(:, 7);
    I           = X(:, 8);
    q           = X(:, 9);
    r           = X(:, 10);
    tax_lab     = X(:, 11);
    tax_r       = X(:, 12);
    transf      = X(:, 13);
    
    Kp          = min(Kp, glob.kmax);
    Bp          = min(Bp, Kp .* glob.bmax);
    
    %______________________________________________________________________
    % Record next period states
    if t < T
        Kt(:, t + 1)   = Kp;
        mut(:, t + 1)  = mup;
        Bt(:, t + 1)   = Bp;
        ratt(:, t + 1) = ratp;
    end
    %______________________________________________________________________
    % Record period t objects
    c_wt(:, t)      = c_w;
    c_bt(:, t)      = c_b;
    Yt(:, t)        = Y;
    Lt(:, t)        = L;
    It(:, t)        = I;
    qt(:, t)        = q;
    rt(:, t)        = r;
    tax_labt(:, t)  = tax_lab;
    tax_rt(:, t)    = tax_r;
    transft(:, t)   = transf;
    if strcmp(options.print, 'Y')
        fprintf('GIRF simulation, period \t%i\n', t);
    end
end

%% Pack
girf.Kt         = Kt;
girf.Bt         = Bt;
girf.mut        = mut;
girf.Zt         = Zt;
girf.c_wt       = c_wt;
girf.c_bt       = c_bt;
girf.Yt         = Yt;
girf.Lt         = Lt;
girf.It         = It;
girf.qt         = qt;
girf.rt         = rt;
girf.tax_labt   = tax_labt;
girf.tax_rt     = tax_rt;
girf.transft    = transft;
girf.X0         = X0;
end