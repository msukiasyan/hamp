function [ys,check] = ngm_steadystate(ys,exo)

% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values fpr the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impos restriction on parameters)

global M_ 

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ 'param.' paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% LOAD PARAMETERS
A = exo(1);
beta = param.beta;
alpha = param.alpha;
delta = param.delta;
gamma = param.gamma;
b_ss = param.b_ss;
chi = param.chi;
nu = param.nu;

R = 1 / beta;
k_l = ((R - (1 - delta)) / (alpha * A)) ^ (1 / (alpha - 1));
W = ((1 - alpha) * A * (k_l) ^ (alpha));

l = fzero(@(x) chi * x^nu - W * ((1 - beta) * b_ss + W * x)^(-gamma), 1);
k = k_l * l;
c_w = (1 - beta) * b_ss + W * l;
c_b = A * k ^ alpha * l ^ (1 - alpha) + ( - delta) * k - c_w;
b = b_ss;


%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = param.' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
