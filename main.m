clc; clear all;
dbstop if error
addpath(genpath(pwd));

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
options.AR1         = 'N';          % Is the shock AR(1) (Y) or only a Markov chain approx (N)
options.disaster    = 'Y';          % Are there diaster shocks? Works only when options.AR1 == 'N'
options.lev_cap     = 'N';          % Is there leverage cap? Defined as d' <= chi * q * k'
options.FB          = 'N';          % If 'Y' then solve_vf_FR computes the first best

% Generalized IRF options
options.Ngirf       = 50000;        % Number of paths to simulate
options.Tgirf       = 50;           % Number of periods to simulate

%% Statespace parameters
glob.n          = [17, 18, 3];      % Number of nodes in each dimension
glob.nf         = [80, 80, 3];      % Number of nodes in each dimension (fine grids)
glob.curv       = [0.2, 2.0];       % Curvature for k (1 is no curvature)
glob.spliorder  = [3, 3, 1];        % Order of splines
glob.kmin       = 1.0;              % Lower bound on capital
glob.kmax       = 40.0;             % Upper bound on capital
glob.bmin       = 0.0;             % Lower bound on deposits
glob.bmax       = 1.00;             % Upper bound on deposits
glob.pzlb       = 0.005;            % Lower bound on probability of z
glob.Ne1        = 10;               % # of approx nodes of AR(1) iid shock in Expectation
glob.plb        = 0.001;            % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2        = 200;              % # of approx nodes of AR(1) iid shock in Approx of Q

%% Full Ramsey statespace parameters
glob.n_FR           = [13, 12, 5, 3];      % Number of nodes in each dimension
glob.n_FR_old       = [13, 12, 5, 3];      % Number of nodes in each dimension
glob.nf_FR          = [20, 80, 20, 3];      % Number of nodes in each dimension (fine grids)
glob.curv_FR        = [0.2, 2.0, 1.2];       % Curvature for k (1 is no curvature)
glob.spliorder_FR   = [3, 3, 3, 1];        % Order of splines

glob.kmin_FR        = 1.0;              % Lower bound on capital
glob.kmax_FR        = 30.0;%25.0             % Upper bound on capital
glob.bmin_FR        = 0.0;             % Lower bound on deposits
glob.bmax_FR        = 0.95;             % Upper bound on deposits
glob.mumin_FR       = -0.05;             % Upper bound on deposits
glob.mumax_FR       = 0.020;%0.004             % Upper bound on deposits
glob.pzlb_FR        = 0.005;            % Lower bound on probability of z
glob.Ne1_FR         = 10;               % # of approx nodes of AR(1) iid shock in Expectation
glob.plb_FR         = 0.001;            % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2_FR         = 200;              % # of approx nodes of AR(1) iid shock in Approx of Q

%% Parameters
glob.beta_w         = 0.97;         % Discount rate of workers
glob.beta_b         = 0.97;         % Discount rate of bankers
glob.beta           = glob.beta_w;  % Uniform discount rate; important for the planner's problem
glob.rhoz           = 0.70;         % AR factor for aggregate productivity
glob.sige           = 0.12;         % Std deviation of the innovation to productivity
glob.dis_prob       = 0.008;        % Conditional prob of entering in disaster state
glob.dis_prod       = 0.70;         % Productivity in disasters
glob.dis_qual       = 0.99; %0.92;        % Productivity in disasters
glob.dis_pers       = 0.66;         % Conditional prob of staying in the disaster state
glob.xi             = 0.80;         % Elasticity parameter of the capital technology
glob.A              = 1 / glob.xi;  % Scale parameter of the capital technology

param.gamma         = 2.0;          % Risk aversion
glob.lambda         = 0.1 ^ param.gamma;          % Pareto weight on experts

param.nu            = 1;%150          % Inverse Frisch
param.nufinal       = 1;
param.Nseq          = 15;
param.nucurv        = 1.5;

param.chi           = 1.0;          % Labor disutility coeff
param.alpha         = 0.36;         % Capital share
param.delta         = 0.025;        % Capital depreciation
param.chi_lev       = 0.8;%1.0178;         % Leverage cap; defined as d' <= chi * q * d'

%% Setup
tic;
[param, glob]       = setup(param, glob, options);
toc;

%% Setup government policy
glob.transf         = 0.0;
glob.bound_pol      = Inf;
glob.dep_tax        = 0.00;

%% Solve for the equilibrium
tic;
eqbm                = solve_eqbm(param, glob, options);
toc;

mean_K              = sum(eqbm.dist .* glob.sf(:, 1));

% Plot
plot_eqbm(eqbm, param, glob, options, mean_K, glob.zgrid(1));

% Compute welfare
[c_v_w, c_v_b]      = compute_value_function(eqbm, param, glob, options);
v_w                 = glob.Phi * c_v_w;
v_b                 = glob.Phi * c_v_b;
v_w_arr             = reshape(v_w, glob.Nk, glob.Nb, glob.Nz);
v_b_arr             = reshape(v_b, glob.Nk, glob.Nb, glob.Nz);

%% Solve simple planner problem
opt_pol_simple          = fmincon(@(x) -worker_welfare_pol(x(1), x(2), [mean_K, 0.8, 0.90], param, glob, options), ...
    [0.85, 0.1], [], [], [], [], [0.6; 0.01], [0.98; 0.15], [], optimoptions('fsolve', 'Display', 'iter'));

% Solve for the equilibrium with policy
glob.target_lev     = 0.8;%opt_pol_simple(1);
glob.slope_pol      = 0.1372;%opt_pol_simple(2);
glob.dep_tax        = 0.00;
glob.transf         = glob.slope_pol * (glob.s(:, 2) >= glob.target_lev) .* (glob.s(:, 3) < 2) .* abs(glob.s(:, 2) - glob.target_lev);
glob.bound_pol      = Inf;
options.save_eqbm   = 'N';
eqbm_pol            = solve_eqbm(param, glob, options);

plot_eqbm_pol(eqbm, eqbm_pol, param, glob, options, mean_K, glob.zgrid(1));

% Compute welfare
[c_v_w_pol, c_v_b_pol]      = compute_value_function(eqbm_pol, param, glob, options);
v_w_pol             = glob.Phi * c_v_w_pol;
v_b_pol             = glob.Phi * c_v_b_pol;
v_w_pol_arr         = reshape(v_w_pol, glob.Nk, glob.Nb, glob.Nz);
v_b_pol_arr         = reshape(v_b_pol, glob.Nk, glob.Nb, glob.Nz);

c_w_arr             = reshape(eqbm.c_w, glob.Nk, glob.Nb, glob.Nz);
c_b_arr             = reshape(eqbm.c_b, glob.Nk, glob.Nb, glob.Nz);
c_w_pol_arr         = reshape(eqbm_pol.c_w, glob.Nk, glob.Nb, glob.Nz);
c_b_pol_arr         = reshape(eqbm_pol.c_b, glob.Nk, glob.Nb, glob.Nz);
% figure;
% subplot(1, 2, 1);
% surf(glob.kgrid, glob.bgrid, max(c_w_pol_arr(:, :, 2)' - c_w_arr(:, :, 2)', -Inf));
% subplot(1, 2, 2);
% surf(glob.kgrid, glob.bgrid, min(c_b_pol_arr(:, :, 2)' - c_b_arr(:, :, 2)', Inf));

[ceq_w, ceq_b]      = cons_equiv(v_w_pol, v_b_pol, eqbm, param, glob, options);
ceq_w_arr = reshape(ceq_w * 100 - 100, glob.Nk, glob.Nb, glob.Nz);
ceq_b_arr = reshape(ceq_b * 100 - 100, glob.Nk, glob.Nb, glob.Nz);
figure;
surf(glob.kgrid, glob.bgrid, ceq_w_arr(:, :, 2)');
figure;
surf(glob.kgrid, glob.bgrid, ceq_b_arr(:, :, 2)');

plot_welfare_simple(ceq_w, ceq_b, mean_K, 0.90, param, glob, options);

% Stats
ratio_CDF           = arrayfun(@(x) sum(eqbm.dist .* (glob.sf(:, 2) <= x)), glob.bgridf);
top1perc            = interp1(ratio_CDF, glob.bgridf, 0.99);
ratio_pol_CDF       = arrayfun(@(x) sum(eqbm_pol.dist .* (glob.sf(:, 2) <= x)), glob.bgridf);
perc_pol            = 100 - 100 * interp1(glob.bgridf, ratio_pol_CDF, top1perc);

% Ex-post bailouts
% expost_bailouts(eqbm, param, glob, options);

%% Solve the full Ramsey problem
% setup with the full Ramsey parameters
tic;
[param, glob_FR]        = setup_FR(param, glob, options);
toc;

options.FB              = 'N';
sol                     = solve_vf_FR(eqbm, glob, glob_FR_old, param, glob_FR, options);

% Used to iteratively update the grids
glob_FR_old             = glob_FR;                                          

%% Solve the first best
options.FB              = 'Y';
sol_FB                  = solve_vf_FR(eqbm, glob, glob_FR_old, param, glob_FR, options);

%% Simulate
sim                     = simulate_model([14, 6, 1], 1, 50000, eqbm, param, glob, options);
sim_pol                 = simulate_model([14, 6, 1], 1, 50000, eqbm_pol, param, glob, options);

%% Recessions
reces                   = find_recessions(sim, param, glob, options);
reces_pol               = find_recessions(sim_pol, param, glob, options, reces.reces_times);

plot_recessions(reces, sim, param, glob, options, reces_pol);

%% 