%% Options
options.Nnewt       = 25;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on policy functions
options.print       = 'N';      % Print out c-solution convergence

%% Statespace parameters
glob.n          = [20, 20, 5];  % Number of nodes in each dimension
glob.curv       = [1.0, 1.0];   % Curvature for k (1 is no curvature)
glob.spliorder  = [3, 3, 3];    % Order of splines (always use linear if productivity is discrete (not AR1))
glob.kmin       = 0.001;        % Lower bound on capital
glob.kmax       = 30.0;         % Upper bound on capital
glob.bmin       = 0.001;        % Lower bound on deposits
glob.bmax       = 35.0;         % Upper bound on deposits
glob.pzlb       = 0.005;        % Lower bound on probability of z
glob.Ne1        = 10;           % # of approx nodes of AR(1) iid shock in Expectation
glob.plb        = 0.001;        % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2        = 200;          % # of approx nodes of AR(1) iid shock in Approx of Q

%% Parameters
glob.beta           = 0.98;     % Patience 
glob.rhoz           = 0.95;     % AR factor for aggregate productivity
glob.sige           = 0.05;     % Std deviation of the innovation to productivity

param.gamma         = 2.0;      % Risk aversion
param.nu            = 1.0;      % Inverse Frisch
param.chi           = 1.0;      % Labor disutility coeff
param.alpha         = 0.36;     % Capital share

%% Setup
tic;
[param, glob]       = setup(param, glob, options);
toc;
%% 