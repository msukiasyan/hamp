function [param, glob] = setup(param, glob, options)

%% State space for aggregate productivity
Nz                  = glob.n(2);                                            % Number of nodes for z
% logzub              = norminv(1 - glob.pzlb, 0, glob.sige / sqrt(1 - glob.rhoz^2));     % log upper bound for z (from its uncond dist)
% logzlb              = norminv(glob.pzlb, 0, glob.sige / sqrt(1 - glob.rhoz^2));         % log lower bound for z
% zlb                 = exp(logzlb);                                        % lower bound for z
% zub                 = exp(logzub);                                        % upper bound for z
% zgrid               = nodeunif(Nz, zlb, zub);                             % 1D uniform grid for z
% Use Rouwenhurst
[P, zgrid, Pssz]        = setup_MarkovZ(Nz, glob.sige, glob.rhoz, 1);       % Pssz is stat dist
zgrid0                  = zgrid;

%% State space for endogenous variable K
Nk              = glob.n(1);                                                % Number of nodes for k
curv            = glob.curv(1);
kgrid           = nodeunif(Nk, glob.kmin .^ curv, glob.kmax .^ curv) .^ (1 / curv);     % Adds curvature
kgrid0          = kgrid;                                                    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli', kgrid, 0, glob.spliorder(1)}, ...         % Create fspace object
                         {'spli', zgrid, 0, glob.spliorder(2)});
sgrid           = funnode(fspace);                                          % Get the new grids
s               = gridmake(sgrid);                                          % Stacked grid points
Ns              = size(s,1);                                                % Total number of points

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
kgrid           = sgrid{1};                                                     
zgrid           = sgrid{2};
Nk              = size(kgrid,1); 
Nz              = size(zgrid,1);

%% Compute expectations matrix
Ne              = glob.Ne1;                             
pvec            = nodeunif(Ne, glob.plb, 1 - glob.plb);                     % Make an equi-spaced grid in probabilities
e               = norminv(pvec, 0, glob.sige);                              % Turn to grid on e
w               = normpdf(e, 0, glob.sige);                                 % Invert normal for shocks
w               = w / sum(w);                                               % Compute pdf of shocks
iNe             = ones(Ne, 1);                           
iNs             = ones(Ns, 1);
gfun            = @(z,e) max(min(exp(glob.rhoz * log(z) + e), max(zgrid)), min(zgrid));     % Constrained to lie within nodes
g               = gfun(kron(s(:, 2), iNe), kron(iNs, e));                    % Compute gfun for all combinations of z and e
Phi             = funbas(fspace, [kron(s(:, 1), iNe), g]);                % For each possible state (s',z') tomorrow
Ikronw          = kron(eye(Ns), w');                                        % Prob matrix
glob.Emat       = Ikronw * Phi;                                             % Expectation of basis over all Ne points
glob.Ematsp     = sparse(glob.Emat);                                        % Create a sparse version

%% Create one time only basis matrices
glob.Phi_Z      = splibas(zgrid0, 0, glob.spliorder(2), s(:, 2));           % Used in Newton computing expected values
Phi_K           = splibas(kgrid0, 0, glob.spliorder(1), s(:, 1));
glob.Phi        = dprod(glob.Phi_Z, Phi_K);                   % Used in Newton updating of c
glob.Phisp      = sparse(glob.Phi);                                         % Create a sparse version (doesn't seem too helpful)
[glob.Phiu, glob.Phil]  = lu(glob.Phisp);
glob.Phiinv     = inv(glob.Phi);
glob.basiscast  = glob.Phiinv * glob.Emat * glob.Phiinv;

%% Set equation solver options
glob.opt        = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-15);

%% Declare additional global variables
glob.kgrid0     = kgrid0;
glob.kgrid      = kgrid;
glob.zgrid0     = zgrid0;
glob.zgrid      = zgrid;
glob.P          = P;
glob.Pssz       = Pssz;
glob.Nz         = Nz;
glob.Nk         = Nk;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;
%__________________________________________________________________________
end
