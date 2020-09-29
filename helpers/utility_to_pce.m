% Convert utility to permanent constant consumption equivalent
function out = utility_to_pce(u, param, glob, options)
    out         = ((1 - glob.beta) * (1 - param.gamma) * u) .^ (1 / (1 - param.gamma));
end