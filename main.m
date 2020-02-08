clc; clear all;
dbstop if error

%% Options
options.Nnewt       = 50;           % Maximum number of Newton steps
options.Nbackw      = 50;           % Maximum number of backward iteration steps
options.tolc        = 1e-2;         % Tolerance on policy functions
options.print       = 'Y';          % Print out convergence information
options.guess       = 'saved';      % Where to get the initial guess ('file' | 'saved' | other (attempts generic guesses))
options.optim       = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'trust-region-dogleg', ...
    'SpecifyObjectiveGradient',true, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-30, ...
    'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIterations', 2000);

%% Statespace parameters
glob.n          = [12, 12, 3];      % Number of nodes in each dimension
glob.curv       = [0.3, 2.0];       % Curvature for k (1 is no curvature)
glob.spliorder  = [3, 3, 1];        % Order of splines
glob.kmin       = 2.0;              % Lower bound on capital
glob.kmax       = 40.0;             % Upper bound on capital
glob.bmin       = 0.05;             % Lower bound on deposits
glob.bmax       = 0.95;             % Upper bound on deposits
glob.pzlb       = 0.005;            % Lower bound on probability of z
glob.Ne1        = 10;               % # of approx nodes of AR(1) iid shock in Expectation
glob.plb        = 0.001;            % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2        = 200;              % # of approx nodes of AR(1) iid shock in Approx of Q

%% Parameters
glob.beta_w         = 0.98;         % Discount rate of workers
glob.beta_b         = 0.97;         % Discount rate of bankers
glob.rhoz           = 0.70;         % AR factor for aggregate productivity
glob.sige           = 0.08;         % Std deviation of the innovation to productivity
glob.xi             = 0.50;         % Elasticity parameter of the capital technology
glob.A              = 1 / glob.xi;  % Scale parameter of the capital technology

param.gamma         = 2.0;          % Risk aversion
param.nu            = 1.0;          % Inverse Frisch
param.chi           = 1.0;          % Labor disutility coeff
param.alpha         = 0.36;         % Capital share
param.delta         = 0.025;        % Capital depreciation

%% Setup
tic;
[param, glob]       = setup(param, glob, options);
toc;

%% Solve for the equilibrium
tic;
eqbm                = solve_eqbm(param, glob, options);
toc;

%% 