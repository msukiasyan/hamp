clc; clear all;
dbstop if error

%% Options
options.Nnewt       = 50;           % Maximum number of Newton steps
options.Nbackw      = 150;          % Maximum number of backward iteration steps
options.itermaxL    = 10000;
options.tolL        = 1e-10;
options.tolc        = 1e-5;         % Tolerance on policy functions
options.print       = 'Y';          % Print out convergence information
options.guess       = 'saved';      % Where to get the initial guess ('file' | 'saved' | other (attempts generic guesses))
options.save_eqbm   = 'Y';
options.guess_FR    = 'saved';      % Where to get the initial guess ('saved' | other (attempts generic guesses))
options.save_eqbm_FR    = 'Y';
options.optim       = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'trust-region-dogleg', ...
    'SpecifyObjectiveGradient',true, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-30, ...
    'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIterations', 2000);
options.burnin      = 1000;          % Number of first periods to ignore in the simulation

% Model options
options.GHH         = 'N';

%% Statespace parameters
glob.n          = [17, 15, 3];      % Number of nodes in each dimension
glob.nf         = [30, 30, 10];      % Number of nodes in each dimension (fine grids)
glob.curv       = [0.2, 2.0];       % Curvature for k (1 is no curvature)
glob.spliorder  = [3, 3, 1];        % Order of splines
glob.kmin       = 1.0;              % Lower bound on capital
glob.kmax       = 60.0;             % Upper bound on capital
glob.bmin       = 0.0;             % Lower bound on deposits
glob.bmax       = 1.00;             % Upper bound on deposits
glob.pzlb       = 0.005;            % Lower bound on probability of z
glob.Ne1        = 10;               % # of approx nodes of AR(1) iid shock in Expectation
glob.plb        = 0.001;            % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2        = 200;              % # of approx nodes of AR(1) iid shock in Approx of Q

%% Full Ramsey statespace parameters
glob.n_FR           = [13, 12, 4, 3];      % Number of nodes in each dimension
glob.n_FR_old       = [13, 12, 4, 3];      % Number of nodes in each dimension
glob.nf_FR          = [7, 7, 3, 3];      % Number of nodes in each dimension (fine grids)
glob.curv_FR        = [0.3, 5.5, 1.2];       % Curvature for k (1 is no curvature)
glob.spliorder_FR   = [2, 2, 2, 1];        % Order of splines

glob.kmin_FR        = 1.0;              % Lower bound on capital
glob.kmax_FR        = 25.0;             % Upper bound on capital
glob.bmin_FR        = 0.0;             % Lower bound on deposits
glob.bmax_FR        = 0.70;             % Upper bound on deposits
glob.mumin_FR       = -0.000;             % Upper bound on deposits
glob.mumax_FR       = 0.004;             % Upper bound on deposits
glob.pzlb_FR        = 0.005;            % Lower bound on probability of z
glob.Ne1_FR         = 10;               % # of approx nodes of AR(1) iid shock in Expectation
glob.plb_FR         = 0.001;            % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2_FR         = 200;              % # of approx nodes of AR(1) iid shock in Approx of Q

%% Parameters
glob.beta_w         = 0.97;         % Discount rate of workers
glob.beta_b         = 0.97;         % Discount rate of bankers
glob.beta           = glob.beta_w;  % Uniform discount rate; important for the planner's problem
glob.rhoz           = 0.70;         % AR factor for aggregate productivity
glob.sige           = 0.08;         % Std deviation of the innovation to productivity
glob.xi             = 0.60;         % Elasticity parameter of the capital technology
glob.A              = 1 / glob.xi;  % Scale parameter of the capital technology
glob.lambda         = 0.5;          % Pareto weight on experts

param.gamma         = 2.0;          % Risk aversion
param.nu            = 170.0;          % Inverse Frisch
param.chi           = 1.0;          % Labor disutility coeff
param.alpha         = 0.36;         % Capital share
param.delta         = 0.025;        % Capital depreciation

%% Setup
tic;
[param, glob]       = setup(param, glob, options);
toc;

%% Setup government policy
glob.transf         = 0.0 * (glob.s(:, 2) > 0.5) .* abs(glob.s(:, 2) - 0.5); %0.3
glob.bound_pol      = Inf;

%% Solve for the equilibrium
tic;
eqbm                = solve_eqbm(param, glob, options);
toc;
% Compute welfare
[c_v_w, c_v_b]      = compute_value_function(eqbm, param, glob, options);
v_w                 = glob.Phi * c_v_w;
v_b                 = glob.Phi * c_v_b;
v_w_arr             = reshape(v_w, glob.Nk, glob.Nb, glob.Nz);
v_b_arr             = reshape(v_b, glob.Nk, glob.Nb, glob.Nz);
Ev_w                = sum(eqbm.dist .* v_w);
Ev_b                = sum(eqbm.dist .* v_b);

% Solve for the equilibrium with policy
glob.transf         = 0.2 * (glob.s(:, 2) > 0.7) .* abs(glob.s(:, 2) - 0.7);
glob.bound_pol      = Inf;
options.save_eqbm   = 'N';
eqbm_pol            = solve_eqbm(param, glob, options);
% Compute welfare
[c_v_w_pol, c_v_b_pol]      = compute_value_function(eqbm_pol, param, glob, options);
v_w_pol             = glob.Phi * c_v_w_pol;
v_b_pol             = glob.Phi * c_v_b_pol;
v_w_pol_arr         = reshape(v_w_pol, glob.Nk, glob.Nb, glob.Nz);
v_b_pol_arr         = reshape(v_b_pol, glob.Nk, glob.Nb, glob.Nz);
Ev_w_pol            = sum(eqbm.dist .* v_w_pol);
Ev_b_pol            = sum(eqbm.dist .* v_b_pol);

figure;
subplot(1, 2, 1);
surf(glob.kgrid, glob.bgrid, max(v_w_pol_arr(:, :, 3)' - v_w_arr(:, :, 3)', -Inf));
subplot(1, 2, 2);
surf(glob.kgrid, glob.bgrid, min(v_b_pol_arr(:, :, 3)' - v_b_arr(:, :, 3)', Inf));

% Ex-post bailouts
% expost_bailouts(eqbm, param, glob, options);

%% Solve the full Ramsey problem
% setup with the full Ramsey parameters
tic;
[param, glob_FR]       = setup_FR(param, glob, options);
toc;

sol                     = solve_FR(eqbm, glob, glob_FR_old, param, glob_FR, options);

%% Simulate
% sim                 = simulate_model([20, 0.7, 1], 1, 50000, eqbm, param, glob, options);
% 
% %% Recessions
% reces               = find_recessions(sim, param, glob, options);

%% 