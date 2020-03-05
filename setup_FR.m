function [param, glob] = setup_FR(param, glob, options)

%% State space for aggregate productivity
Nz                  = glob.n_FR(4);                                            % Number of nodes for z
% logzub              = norminv(1 - glob.pzlb, 0, glob.sige / sqrt(1 - glob.rhoz^2));     % log upper bound for z (from its uncond dist)
% logzlb              = norminv(glob.pzlb, 0, glob.sige / sqrt(1 - glob.rhoz^2));         % log lower bound for z
% zlb                 = exp(logzlb);                                        % lower bound for z
% zub                 = exp(logzub);                                        % upper bound for z
% zgrid               = nodeunif(Nz, zlb, zub);                             % 1D uniform grid for z
% Use Rouwenhurst
[P, zgrid, Pssz]        = setup_MarkovZ(Nz, glob.sige, glob.rhoz, 1);       % Pssz is stat dist
zgrid0                  = zgrid;

%% State space for endogenous variable K
Nk              = glob.n_FR(1);                                                % Number of nodes for k
curv            = glob.curv_FR(1);
kgrid           = nodeunif(Nk, glob.kmin_FR .^ curv, glob.kmax_FR .^ curv) .^ (1 / curv);     % Adds curvature
kgrid0          = kgrid;                                                    % Save for computing basis matrices

%% State space for endogenous variable B
Nb              = glob.n_FR(2);                                                % Number of nodes for k
curv            = glob.curv_FR(2);
bgrid           = nodeunif(Nb, glob.bmin_FR .^ curv, glob.bmax_FR .^ curv) .^ (1 / curv);    % Adds curvature
bgrid0          = bgrid;                                                    % Save for computing basis matrices

%% State space for endogenous variable mu
Nmu             = glob.n_FR(3);                                                % Number of nodes for k
curv            = glob.curv_FR(3);
mugrid          = nodeunif(Nmu, glob.mumin_FR .^ curv, glob.mumax_FR .^ curv) .^ (1 / curv);    % Adds curvature
mugrid0         = mugrid;                                                    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli', kgrid, 0, glob.spliorder_FR(1)}, ...         % Create fspace object
                         {'spli', bgrid, 0, glob.spliorder_FR(2)}, ...
                         {'spli', mugrid, 0, glob.spliorder_FR(3)}, ...
                         {'spli', zgrid, 0, glob.spliorder_FR(4)});
sgrid           = funnode(fspace);                                          % Get the new grids
s               = gridmake(sgrid);                                          % Stacked grid points
Ns              = size(s, 1);                                                % Total number of points

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
kgrid           = sgrid{1};                                                     
bgrid           = sgrid{2};
mugrid          = sgrid{3};
zgrid           = sgrid{4};
Nk              = size(kgrid,1); 
Nb              = size(bgrid,1); 
Nmu             = size(mugrid,1); 
Nz              = size(zgrid,1);

%% Compute expectations matrix
Ne              = glob.Ne1_FR;                             
pvec            = nodeunif(Ne, glob.plb_FR, 1 - glob.plb_FR);                     % Make an equi-spaced grid in probabilities
e               = norminv(pvec, 0, glob.sige);                              % Turn to grid on e
w               = normpdf(e, 0, glob.sige);                                 % Invert normal for shocks
w               = w / sum(w);                                               % Compute pdf of shocks
iNe             = ones(Ne, 1);                           
iNs             = ones(Ns, 1);
gfun            = @(z,e) max(min(exp(glob.rhoz * log(z) + e), max(zgrid)), min(zgrid));     % Constrained to lie within nodes
g               = gfun(kron(s(:, 4), iNe),kron(iNs, e));                    % Compute gfun for all combinations of z and e
Phi             = funbas(fspace, [kron(s(:, 1:3), iNe), g]);                % For each possible state (s',z') tomorrow
Ikronw          = kron(eye(Ns), w');                                        % Prob matrix
glob.Emat       = Ikronw * Phi;                                             % Expectation of basis over all Ne points
glob.Ematsp     = sparse(glob.Emat);                                        % Create a sparse version

%% Construct fine grid for histogram
% kgridf          = nodeunif(glob.nf(1), glob.kmin .^ glob.curv(1), glob.kmax .^ glob.curv(1)) .^ (1/glob.curv(1));
% bgridf          = nodeunif(glob.nf(2), glob.bmin .^ glob.curv(2), glob.bmax .^ glob.curv(2)) .^ (1/glob.curv(2));
kgridf          = kgrid;
bgridf          = bgrid;
mugridf         = mugrid;
Nkf             = size(kgridf, 1);
Nbf             = size(bgridf, 1);
Nmuf            = size(mugridf, 1);

% zgridf          = nodeunif(glob.nf(3), min(zgrid), max(zgrid));      
zgridf          = zgrid;
Nzf             = size(zgridf, 1);
sf              = gridmake(kgridf, bgridf, mugridf, zgridf);
Nsf             = size(sf, 1);

glob.kgridf     = kgridf;
glob.bgridf     = bgridf;
glob.mugridf    = mugridf;
glob.zgridf     = zgridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QZ matrix for approximation of stationary distribution
Ne              = glob.Ne2_FR;
pvec            = nodeunif(Ne, glob.plb_FR, 1 - glob.plb_FR);         % Make an equi-spaced grid in probabilities
e               = norminv(pvec, 0, glob.sige);                % Invert normal for shocks
w               = normpdf(e, 0, glob.sige);                   % Compute pdf of shocks
w               = w / sum(w);                                 % Normalise
fspaceZ         = fundef({'spli', zgridf, 0, 1});              % Linear interpolant
QZ              = zeros(Nsf, Nzf);
P               = zeros(Nzf, Nzf);                           % P constructed so can compute steady state Psszf and compare to Pssz
for i = 1:Ne
    g           = gfun(sf(:, 4), e(i));                       % z' for all z and realization e(i)
    QZi         = funbas(fspaceZ, g);                        % Basis for all such z'
    QZ          = QZ + w(i) * QZi;                            % Weight by prob w(i) and add
    P           = P  + w(i) * funbas(fspaceZ, gfun(zgridf, e(i))); % Compute the same thing as above but only for z-states
end
glob.QZ         = QZ;
% For plotting
Psszf           = P ^ 1000;
Psszf           = Psszf(1,:)';
glob.Psszf      = Psszf;                                    % Stat dist for z

%% Create one time only basis matrices
glob.Phi_Z      = splibas(zgrid0, 0, glob.spliorder_FR(4), s(:, 4));           % Used in Newton computing expected values
glob.Phi_Zf     = splibas(zgrid0, 0, glob.spliorder_FR(4), sf(:,4)); 
Phi_K           = splibas(kgrid0, 0, glob.spliorder_FR(1), s(:, 1));
Phi_B           = splibas(bgrid0, 0, glob.spliorder_FR(2), s(:, 2));
Phi_mu          = splibas(mugrid0, 0, glob.spliorder_FR(3), s(:, 3));
glob.Phi        = dprod(glob.Phi_Z, dprod(Phi_mu, dprod(Phi_B, Phi_K)));                   % Used in Newton updating of c
glob.Phisp      = sparse(glob.Phi);                                         % Create a sparse version (doesn't seem too helpful)
[glob.Phiu, glob.Phil]  = lu(glob.Phisp);
glob.Phiinv     = inv(glob.Phi);
glob.basiscast  = glob.Phiinv * glob.Emat * glob.Phiinv;

%% Declare additional global variables
glob.kgrid0     = kgrid0;
glob.kgrid      = kgrid;
glob.bgrid0     = bgrid0;
glob.bgrid      = bgrid;
glob.mugrid0    = mugrid0;
glob.mugrid     = mugrid;
glob.zgrid0     = zgrid0;
glob.zgrid      = zgrid;
glob.P          = P;
glob.Pssz       = Pssz;
glob.Nz         = Nz;
glob.Nk         = Nk;
glob.Nb         = Nb;
glob.Nmu        = Nmu;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;
%__________________________________________________________________________
end
