function [P, zgrid, Pssz]   = disaster_Markov(P0, zgrid0, Pssz0, param, glob, options)
%% Unpack
dis_prob        = glob.dis_prob;
dis_prod        = glob.dis_prod;
dis_pers        = glob.dis_pers;
Nz0             = size(zgrid0, 1);

%% Modify the transition matrix and grid to add the disaster state
% zgrid           = [zgrid0(1) - 1e-9; zgrid0];
zgrid           = [dis_prod; zgrid0];
P               = [dis_pers, (1 - dis_pers) * Pssz0';
                    dis_prob * ones(Nz0, 1), (1 - dis_prob) * P0];
                
Pssz            = P^10000;
Pssz            = Pssz(1, :)';

end