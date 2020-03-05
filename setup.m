function [param, glob] = setup(param, glob, options)

%% State space for aggregate productivity
Nz                  = glob.n(3);                                            % Number of nodes for z
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

%% State space for endogenous variable B
Nb              = glob.n(2);                                                % Number of nodes for k
curv            = glob.curv(2);
bgrid           = nodeunif(Nb, glob.bmin .^ curv, glob.bmax .^ curv) .^ (1 / curv);    % Adds curvature
bgrid0          = bgrid;                                                    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli', kgrid, 0, glob.spliorder(1)}, ...         % Create fspace object
                         {'spli', bgrid, 0, glob.spliorder(2)}, ...
                         {'spli', zgrid, 0, glob.spliorder(3)});
sgrid           = funnode(fspace);                                          % Get the new grids
s               = gridmake(sgrid);                                          % Stacked grid points
Ns              = size(s,1);                                                % Total number of points

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
kgrid           = sgrid{1};                                                     
bgrid           = sgrid{2};
zgrid           = sgrid{3};
Nk              = size(kgrid,1); 
Nb              = size(bgrid,1); 
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
g               = gfun(kron(s(:, 3), iNe),kron(iNs, e));                    % Compute gfun for all combinations of z and e
Phi             = funbas(fspace, [kron(s(:, 1:2), iNe), g]);                % For each possible state (s',z') tomorrow
Ikronw          = kron(eye(Ns), w');                                        % Prob matrix
glob.Emat       = Ikronw * Phi;                                             % Expectation of basis over all Ne points
glob.Ematsp     = sparse(glob.Emat);                                        % Create a sparse version

%% Construct fine grid for histogram
kgridf2         = nodeunif(glob.nf(1), glob.kmin .^ glob.curv(1), glob.kmax .^ glob.curv(1)) .^ (1/glob.curv(1));
bgridf2         = nodeunif(glob.nf(2), glob.bmin .^ glob.curv(2), glob.bmax .^ glob.curv(2)) .^ (1/glob.curv(2));
kgridf          = kgrid;
bgridf          = bgrid;
Nkf             = size(kgridf, 1);
Nbf             = size(bgridf, 1);

zgridf2         = nodeunif(glob.nf(3), min(zgrid), max(zgrid));      
zgridf          = zgrid;
Nzf             = size(zgridf, 1);
sf              = gridmake(kgridf, bgridf, zgridf);
Nsf             = size(sf, 1);

glob.kgridf     = kgridf;
glob.bgridf     = bgridf;
glob.zgridf     = zgridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QZ matrix for approximation of stationary distribution
Ne              = glob.Ne2;
pvec            = nodeunif(Ne, glob.plb, 1 - glob.plb);         % Make an equi-spaced grid in probabilities
e               = norminv(pvec, 0, glob.sige);                % Invert normal for shocks
w               = normpdf(e, 0, glob.sige);                   % Compute pdf of shocks
w               = w / sum(w);                                 % Normalise
fspaceZ         = fundef({'spli', zgridf, 0, 1});              % Linear interpolant
QZ              = zeros(Nsf, Nzf);
P               = zeros(Nzf, Nzf);                           % P constructed so can compute steady state Psszf and compare to Pssz
for i = 1:Ne
    g           = gfun(sf(:, 3), e(i));                       % z' for all z and realization e(i)
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
glob.Phi_Z      = splibas(zgrid0, 0, glob.spliorder(3), s(:, 3));           % Used in Newton computing expected values
glob.Phi_Zf     = splibas(zgrid0, 0, glob.spliorder(3), sf(:,3)); 
Phi_K           = splibas(kgrid0, 0, glob.spliorder(1), s(:, 1));
Phi_B           = splibas(bgrid0, 0, glob.spliorder(2), s(:, 2));
glob.Phi        = dprod(glob.Phi_Z, dprod(Phi_B, Phi_K));                   % Used in Newton updating of c
glob.Phisp      = sparse(glob.Phi);                                         % Create a sparse version (doesn't seem too helpful)
[glob.Phiu, glob.Phil]  = lu(glob.Phisp);
glob.Phiinv     = inv(glob.Phi);
glob.basiscast  = glob.Phiinv * glob.Emat * glob.Phiinv;

%% Declare additional global variables
glob.kgrid0     = kgrid0;
glob.kgrid      = kgrid;
glob.bgrid0     = bgrid0;
glob.bgrid      = bgrid;
glob.zgrid0     = zgrid0;
glob.zgrid      = zgrid;
glob.P          = P;
glob.Pssz       = Pssz;
glob.Nz         = Nz;
glob.Nk         = Nk;
glob.Nb         = Nb;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;
glob.kgridf2    = kgridf2;
glob.bgridf2    = bgridf2;
glob.zgridf2    = zgridf2;
glob.sf2        = gridmake(kgridf2, bgridf2, zgridf2);
%__________________________________________________________________________
end
