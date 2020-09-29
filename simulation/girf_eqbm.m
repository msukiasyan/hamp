function girf = girf_eqbm(s0, eq, param, glob, options)
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
ratt(:, 1)  = st(2) / st(1);

% Policy interpolants
cKp             = Phil \ (Phiu \ eq.Kp);
cBp             = Phil \ (Phiu \ eq.Bp);
cc_w            = Phil \ (Phiu \ eq.c_w);
cc_b            = Phil \ (Phiu \ eq.c_b);
cL              = Phil \ (Phiu \ eq.L);
cY              = Phil \ (Phiu \ eq.Y);
cI              = Phil \ (Phiu \ eq.I);
cq              = Phil \ (Phiu \ eq.q);
cr              = Phil \ (Phiu \ eq.r);

% Draw paths from the Markov chain
rng(999);
Zi          = find(zgrid == s0(3));
simcnt      = zeros(glob.Nz, 1);
simcnt(Zi)  = Nsim;
% simcnt(Zi)                  = Nsim / 2;
% simcnt(median(glob.Nz))     = Nsim / 2;
Zi          = simulate(dtmc(P), T - 1, 'X0', simcnt)';
Zt          = zgrid(Zi);

% Compute all variables if the shock was at the ergodic mode
X0          = funeval([cKp, cBp, cc_w, cc_b, cY, cL, cI, cq, cr], fspace, [st(1), st(2) / st(1), median(zgrid)]);
Kp0         = X0(1);
Bp0         = X0(2);
c_w0        = X0(3);
c_b0        = X0(4);
Y0          = X0(5);
L0          = X0(6);
I0          = X0(7);
q0          = X0(8);
r0          = X0(9);

% Simulation
for t = 1:T
    % Current state
    K       = Kt(:, t);
    B       = Bt(:, t);
    Z       = Zt(:, t);
    rat     = ratt(:, t);
    % Compute policies given states
    X       = funeval([cKp, cBp, cc_w, cc_b, cY, cL, cI, cq, cr], fspace, [K, rat, Z]);
    
    Kp          = X(:, 1);
    Bp          = X(:, 2);
    
    ratp        = Bp ./ Kp;
    qual_shock  = (glob.dis_qual * (Zi == 1) + 1 * (Zi ~= 1));
    K           = K .* qual_shock;
    
    c_w         = X(:, 3);
    c_b         = X(:, 4);
    Y           = X(:, 5);
    L           = X(:, 6);
    I           = X(:, 7);
    q           = X(:, 8);
    r           = X(:, 9);
    
    Kp          = min(Kp, glob.kmax);
    Bp          = min(Bp, Kp .* glob.bmax);
    
    %______________________________________________________________________
    % Record next period states
    if t < T
        Kt(:, t + 1)   = Kp;
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
    if strcmp(options.print, 'Y')
        fprintf('GIRF simulation, period %i\n', t);
    end
end

%% Pack
girf.Kt         = Kt;
girf.Bt         = Bt;
girf.Zt         = Zt;
girf.c_wt       = c_wt;
girf.c_bt       = c_bt;
girf.Yt         = Yt;
girf.Lt         = Lt;
girf.It         = It;
girf.qt         = qt;
girf.rt         = rt;

girf.Kp0        = Kp0;
girf.Bp0        = Bp0;
girf.c_w0       = c_w0;
girf.c_b0       = c_b0;
girf.Y0         = Y0;
girf.L0         = L0;
girf.I0         = I0;
girf.q0         = q0;
girf.r0         = r0;
end