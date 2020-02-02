function [Pz,zvec,Pssz] = setup_MarkovZ(Nz,sige,rhoz,A)
% by Simon Mongey
%________________________________________________________________________
% Function
% 1. Using rouwenhorst
sig_uncond      = sige/sqrt(1-rhoz^2); % Uncond variance of log z
mu_uncond       = log(A)-(1/2)*(sige^2)/((1+rhoz)*(1-rhoz)); % Uncond mean of log z
[Pz, logzvec]   = setup_rouwen(rhoz,mu_uncond,sig_uncond,Nz);  % Pz is the transition matrix, logzvec is the grid
Pz              = Pz'; % Pz(i, j) is Prob(y'=z_j | y=z_i)
zvec            = exp(logzvec);
Pssz            = Pz^10000; % Repeat for stat dist
Pssz            = Pssz(1,:)'; % All rows contain the stat dist; take the first one
% 2. Remove points with Pr(z)<1e-5 in the stationary distribution
% Issue: The rouwen procedure creates many points with little prob
cut = 1e-6;       % Chop off points with low productivity
if (cut>1e-8)&&(min(Pssz)<cut)
    go  = 1;
    Nztemp = Nz;
    while go % Increase Nztemp until the stat dist has at least Nz non-small probs
        Nztemp          = Nztemp+1;
        [Pz, logzvec]   = setup_rouwen(rhoz,mu_uncond,sig_uncond,Nztemp);
        Pz              = Pz';
        zvec            = exp(logzvec);
        Pssz            = Pz^10000;
        Pssz            = Pssz(1,:)';
        X               = sum(Pssz>cut);
        if X>=Nz;go=0;end;
    end
    i       = (Pssz>cut); % Keep the non-small prob gridpoints only
    zvec    = zvec(i);
    Pz      = Pz(i,i);
end
Pz          = bsxfun(@rdivide,Pz,sum(Pz,2)); % Normalize so that rows sum to 1
Pz(:,end)   = 1-sum(Pz(:,1:end-1),2);
Pssz        = Pz^10000; % Recompute stat dist
Pssz        = Pssz(1,:)';

