dynare ngm nostrict;

bmin = -15.0;
bmax = 30;
Nb = 50;
bgrid = bmin + (linspace(0, 1, Nb) .^ 0.5) * (bmax - bmin);

kmin = 0.2;
kmax = 1.5;
Nk = 15;
kgrid = kmin + (linspace(0, 1, Nk) .^ 2) * (kmax - kmin);

kgrid = [kgrid, linspace(53, 55, 5)];
Nk  = Nk + 5;

res = cell(Nb * Nk, 1);

fres = zeros(201 * Nb * Nk, 5);

for bi = 1:Nb
    for ki = 1:Nk
        nn = (bi - 1) * Nk + ki;
        set_param_value('b_ss', bgrid(bi));
        disp(['Trying b_ss = ' num2str(bgrid(bi))]);
        evalc('perfect_foresight_setup;');
        oo_.endo_simul(3, 1) = kgrid(ki);
        evalc('steady; perfect_foresight_solver;');
        res{nn} = [k(1:end-1), l(2:end), W(2:end), b(2:end), c_w(2:end), c_b(2:end)];
        fres((nn - 1) * 201 + 1:nn * 201, 1) = k(1:end-1);
        fres((nn - 1) * 201 + 1:nn * 201, 2) = b(2:end);
        fres((nn - 1) * 201 + 1:nn * 201, 3) = c_w(2:end);
        fres((nn - 1) * 201 + 1:nn * 201, 4) = c_b(2:end);
        fres((nn - 1) * 201 + 1:nn * 201, 5) = R(2:end);
    end
end

figure;
subplot(1, 2, 1); scatter3(fres(:, 1), fres(:, 2), fres(:, 3));
title('"Worker"');
xlabel('K');
ylabel('B');
zlabel('C');

subplot(1, 2, 2); scatter3(fres(:, 1), fres(:, 2), fres(:, 4));
title('"Banker"');
xlabel('K');
ylabel('B');
zlabel('C');

save('..\init_guess.mat', 'fres');
