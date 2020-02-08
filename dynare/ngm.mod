var c_w c_b k l R W b;
varexo A;
parameters alpha beta delta gamma b_ss chi nu k0;

alpha = 0.36;
beta = 0.97;
delta = 0.025;
gamma = 1;
b_ss = 1;
chi = 1;
nu = 1;
k0 = 0.2;

model;
c_w^(-gamma) = beta * R(+1) * c_w(+1)^(-gamma);
c_b^(-gamma) = beta * R(+1) * c_b(+1)^(-gamma);

chi * l^nu / c_w^(-gamma) = W;

R = (alpha * A * (k(-1)/l) ^ (alpha - 1) + 1 - delta);
W = ((1 - alpha) * A * (k(-1)/l) ^ (alpha));
c_w + c_b + k = A * k(-1) ^ alpha * l ^ (1 - alpha) + (1 - delta) * k(-1);
b(+1)  = R * (b - c_w + W * l);
end;

initval;
A=1;
k = k0;% 0.5; %0.5*((1-beta*(1-delta))/(beta*alpha*A))^(1/(alpha-1));

end;

endval;
A = 1;
end;
steady;

% histval;
% end;

simul(periods=200);

% rplot k;
% rplot c_w ;
% rplot c_b;
% rplot b ;