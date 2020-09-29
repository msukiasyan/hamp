function plot_planner_welfare(sol, eq, glob_eq, param, glob, options, Kfix, mufix, Zfix)
%% Unpack
c_w         = sol.c_w;
c_b         = sol.c_b;
L           = sol.L;
Y           = sol.Y;
I           = sol.I;
q           = sol.q;
Kp          = sol.Kp;
Bp          = sol.Bp;
r           = sol.r;
cold        = sol.cold;
ecold       = sol.ecold;
dist        = sol.dist;
s           = glob.s;
Phi         = glob.Phi;
Phiu        = glob.Phiu;
Phil        = glob.Phil;

%% Plot welfare
ps          = gridmake([Kfix * 0.9; Kfix; Kfix * 1.1], glob.bgridf, mufix, [Zfix * 0.9; Zfix; Zfix * 1.1]);

Phip        = funbas(glob.fspace, ps, [0, 0, 0, 0]);
Phip_eq     = funbas(glob_eq.fspace, ps(:, [1, 2, 4]), [0, 0, 0]);

[c_v_w_pl, c_v_b_pl]    = compute_value_pol(sol, param, glob, options);
[c_v_w_eq, c_v_b_eq]    = compute_value_function(eq, param, glob_eq, options);

v_w_pl      = Phip * c_v_w_pl;
v_b_pl      = Phip * c_v_b_pl;
v_w_eq      = Phip_eq * c_v_w_eq;
v_b_eq      = Phip_eq * c_v_b_eq;

v_w_pl_arr  = reshape(v_w_pl, 3, glob.Nbf, 3);
v_b_pl_arr  = reshape(v_b_pl, 3, glob.Nbf, 3);
v_w_eq_arr  = reshape(v_w_eq, 3, glob.Nbf, 3);
v_b_eq_arr  = reshape(v_b_eq, 3, glob.Nbf, 3);

figure;
subplot(1, 2, 1);
plot(glob.bgridf, v_w_pl_arr(:, :, 1), 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, v_w_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('D / K')
ylabel('C_w')
legend(['K at -10%'], ['K at mean'], ['K at +10%']);
title('Worker welfare, Z at -10%')
grid on

subplot(1, 2, 2);
plot(glob.bgridf, v_b_pl_arr(:, :, 1), 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, v_b_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('D / K')
ylabel('C_e')
legend(['K at -10%'], ['K at mean'], ['K at +10%']);
title('Expert welfare, Z at -10%')
grid on


end