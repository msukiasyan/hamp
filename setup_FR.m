function [param, glob_FR] = setup_FR(param, glob, options)
glob_FR                 = glob;

%% State space for aggregate productivity
Nz                      = glob.n_FR(4);
switch options.AR1
    case 'Y' 
        % Use Rouwenhurst
        [P, zgrid, Pssz]    = setup_MarkovZ(Nz, glob.sige, glob.rhoz, 1);       % Pssz is stat dist
    case 'N'
        [P, zgrid, Pssz]    = setup_MarkovZ(Nz, glob.sige, glob.rhoz, 1);         
        if options.disaster == 'Y'
            [P, zgrid, Pssz]    = disaster_Markov(P, zgrid, Pssz, param, glob, options);
        end
end
zgrid0                  = zgrid;

%% State space for endogenous variable K
Nk              = glob.n_FR(1);                                                % Number of nodes for k
curv            = glob.curv_FR(1);
kgrid           = nodeunif(Nk, glob.kmin_FR .^ curv, glob.kmax_FR .^ curv) .^ (1 / curv);     % Adds curvature

% add extra points
% kgrid           = sort([kgrid; (kgrid(Nk-4:Nk) + kgrid(Nk-5:Nk-1)) / 2]);

% extragrid       = nodeunif(5, ((glob.kmin_FR + glob.kmax_FR) / 2) .^ 5.0, (glob.kmax_FR - 0.5) .^ 5.0) .^ (1 / 5.0);     % Adds curvature
% kgrid           = sort([kgrid; extragrid]);
% Nk              = Nk + 5;
% glob_FR.n_FR(1) = Nk;

kgrid0          = kgrid;                                                    % Save for computing basis matrices

%% State space for endogenous variable B
Nb              = glob.n_FR(2);                                                % Number of nodes for k
curv            = glob.curv_FR(2);
bgrid           = nodeunif(Nb, glob.bmin_FR .^ curv, glob.bmax_FR .^ curv) .^ (1 / curv);    % Adds curvature
bgrid0          = bgrid;                                                    % Save for computing basis matrices

%% State space for endogenous variable mu
Nmu             = glob.n_FR(3);                                                % Number of nodes for k
curv            = glob.curv_FR(3);
mugrid          = nodeunif(Nmu, 0, (glob.mumax_FR - glob.mumin_FR) .^ curv) .^ (1 / curv) + glob.mumin_FR;    % Adds curvature
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
switch options.AR1
    case 'Y'
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
        glob_FR.Emat    = Ikronw * Phi;                                             % Expectation of basis over all Ne points
    case 'N'
        Phi             = funbas(fspace, s);                                % Allowing only z' on zgrid
        glob_FR.Emat    = kron(P, speye(Nk * Nb * Nmu)) * Phi;   
end

glob_FR.Ematsp     = sparse(glob_FR.Emat);                                        % Create a sparse version

%% Construct fine grid for histogram
kgridf          = nodeunif(glob.nf_FR(1), glob.kmin_FR .^ glob.curv_FR(1), glob.kmax_FR .^ glob.curv_FR(1)) .^ (1/glob.curv_FR(1));
bgridf          = nodeunif(glob.nf_FR(2), glob.bmin_FR .^ glob.curv_FR(2), glob.bmax_FR .^ glob.curv_FR(2)) .^ (1/glob.curv_FR(2));
mugridf         = nodeunif(glob.nf_FR(3), 0, (glob.mumax_FR - glob.mumin_FR) .^ curv) .^ (1 / curv) + glob.mumin_FR;    % Adds curvature
Nkf             = size(kgridf, 1);
Nbf             = size(bgridf, 1);
Nmuf            = size(mugridf, 1);

switch options.AR1
    case 'Y'
        zgridf      = nodeunif(glob.nf_FR(3), min(zgrid), max(zgrid));
    case 'N'
        zgridf      = zgrid;                                                % Ignore glob.nf(3) because have discrete z's
end 

Nzf             = size(zgridf, 1);
sf              = gridmake(kgridf, bgridf, mugridf, zgridf);
Nsf             = size(sf, 1);

%% Compute QZ matrix for approximation of stationary distribution
switch options.AR1
    case 'Y' 
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
        glob_FR.QZ         = QZ;
        % For plotting
        Psszf           = P ^ 1000;
        Psszf           = Psszf(1,:)';
        glob_FR.Psszf      = Psszf;                                    % Stat dist for z
    case 'N'
        Pf              = P;
        glob_FR.QZ      = kron(P, ones(Nkf * Nbf * Nmuf, 1));  
end

%% Create one time only basis matrices
glob_FR.Phi_Z       = splibas(zgrid0, 0, glob.spliorder_FR(4), s(:, 4));           % Used in Newton computing expected values
glob_FR.Phi_Zf      = splibas(zgrid0, 0, glob.spliorder_FR(4), sf(:,4)); 
Phi_K               = splibas(kgrid0, 0, glob.spliorder_FR(1), s(:, 1));
Phi_B               = splibas(bgrid0, 0, glob.spliorder_FR(2), s(:, 2));
Phi_mu              = splibas(mugrid0, 0, glob.spliorder_FR(3), s(:, 3));
glob_FR.Phi         = dprod(glob_FR.Phi_Z, dprod(Phi_mu, dprod(Phi_B, Phi_K)));                   % Used in Newton updating of c
glob_FR.Phisp       = sparse(glob_FR.Phi);                                         % Create a sparse version (doesn't seem too helpful)
[glob_FR.Phiu, glob_FR.Phil]  = lu(glob_FR.Phisp);
glob_FR.Phif        = funbas(fspace, sf);
% glob_FR.Phiinv      = inv(glob_FR.Phi);
% glob_FR.basiscast   = glob_FR.Phiinv * glob_FR.Emat * glob_FR.Phiinv;
% glob_FR.PhiinvEmat  = glob_FR.Phiinv * glob_FR.Emat;

%% Declare additional global variables
glob_FR.kgrid0     = kgrid0;
glob_FR.kgrid      = kgrid;
glob_FR.bgrid0     = bgrid0;
glob_FR.bgrid      = bgrid;
glob_FR.mugrid0    = mugrid0;
glob_FR.mugrid     = mugrid;
glob_FR.zgrid0     = zgrid0;
glob_FR.zgrid      = zgrid;
glob_FR.P          = P;
glob_FR.Pf         = Pf;
glob_FR.Pssz       = Pssz;
glob_FR.Nz         = Nz;
glob_FR.Nk         = Nk;
glob_FR.Nb         = Nb;
glob_FR.Nmu        = Nmu;
glob_FR.fspace     = fspace;
glob_FR.s          = s;
glob_FR.Ns         = Ns;
glob_FR.kgridf     = kgridf;
glob_FR.bgridf     = bgridf;
glob_FR.mugridf    = mugridf;
glob_FR.zgridf     = zgridf;
glob_FR.sf         = sf;
glob_FR.Nsf        = Nsf;
glob_FR.Nkf        = Nkf;
glob_FR.Nbf        = Nbf;
glob_FR.Nmuf       = Nmuf;
glob_FR.Nzf        = Nzf;

%__________________________________________________________________________
end
