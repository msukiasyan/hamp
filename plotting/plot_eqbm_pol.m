function plot_eqbm_pol(eq, eq_pol, param, glob, options, Kfix, Zfix)
%% Unpack
c_w         = eq.c_w;
c_b         = eq.c_b;
L           = eq.L;
Y           = eq.Y;
I           = eq.I;
q           = eq.q;
Kp          = eq.Kp;
Bp          = eq.Bp;
r           = eq.r;
c           = eq.c;
dist        = eq.dist;
lev         = eq.lev;
s           = glob.s;
Phi         = glob.Phi;
Phiu        = glob.Phiu;
Phil        = glob.Phil;

c_w_pol     = eq_pol.c_w;
c_b_pol     = eq_pol.c_b;
L_pol       = eq_pol.L;
Y_pol       = eq_pol.Y;
I_pol       = eq_pol.I;
q_pol       = eq_pol.q;
Kp_pol      = eq_pol.Kp;
Bp_pol      = eq_pol.Bp;
r_pol       = eq_pol.r;
c_pol       = eq_pol.c;
dist_pol    = eq_pol.dist;
lev_pol     = eq_pol.lev;

%% Plot
ps              = gridmake([Kfix * 0.9; Kfix; Kfix * 1.1], glob.bgridf, Zfix);

Phip            = funbas(glob.fspace, ps);

If              = Phip * (Phil \ (Phiu \ I));
If_arr          = reshape(If, 3, glob.Nbf);
If_pol          = Phip * (Phil \ (Phiu \ I_pol));
If_pol_arr      = reshape(If_pol, 3, glob.Nbf);

Lf              = Phip * (Phil \ (Phiu \ L));
Lf_arr          = reshape(Lf, 3, glob.Nbf);
Lf_pol          = Phip * (Phil \ (Phiu \ L_pol));
Lf_pol_arr      = reshape(Lf_pol, 3, glob.Nbf);

rf              = Phip * (Phil \ (Phiu \ r));
rf_arr          = reshape(rf, 3, glob.Nbf);
rf_pol          = Phip * (Phil \ (Phiu \ r_pol));
rf_pol_arr      = reshape(rf_pol, 3, glob.Nbf);

c_wf            = Phip * (Phil \ (Phiu \ c_w));
c_wf_arr        = reshape(c_wf, 3, glob.Nbf);
c_wf_pol        = Phip * (Phil \ (Phiu \ c_w));
c_wf_pol_arr    = reshape(c_wf_pol, 3, glob.Nbf);

c_bf            = Phip * (Phil \ (Phiu \ c_b));
c_bf_arr        = reshape(c_bf, 3, glob.Nbf);
c_bf_pol        = Phip * (Phil \ (Phiu \ c_b));
c_bf_pol_arr    = reshape(c_bf_pol, 3, glob.Nbf);

figure;
subplot(1, 3, 1);
plot(glob.bgridf, If_pol_arr, 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, If_arr, 'LineStyle', '--', 'LineWidth', 2)
xline(glob.target_lev, '--k');
xlabel('$D / K$')
ylabel('$I$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Investment, $Z = ', num2str(Zfix), '$'])
grid on

subplot(1, 3, 2);
plot(glob.bgridf, rf_pol_arr, 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, rf_arr, 'LineStyle', '--', 'LineWidth', 2)
xline(glob.target_lev, '--k');
xlabel('$D / K$')
ylabel('$r$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Risk-free rate, $Z = ', num2str(Zfix), '$'])
grid on

subplot(1, 3, 3);
plot(glob.bgridf, Lf_pol_arr, 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, Lf_arr, 'LineStyle', '--', 'LineWidth', 2)
xline(glob.target_lev, '--k');
xlabel('$D / K$')
ylabel('$L$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Labor, $Z = ', num2str(Zfix), '$'])
grid on

bail_perc           = glob.transf ./ production(glob.s(:, 3), glob.s(:, 1), L_pol, param, glob, options) * 100;
bail_percf          = Phip * (Phil \ (Phiu \ bail_perc));
bail_percf_arr      = reshape(bail_percf, 3, glob.Nbf);

figure;
plot(glob.bgridf, bail_percf_arr, 'LineWidth', 2)
xline(glob.target_lev, '--k');
xlabel('$D / K$')
ylabel('$T / Y$, \%')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Bailouts as percentage of GDP, $Z = ', num2str(Zfix), '$'])
grid on

% Plots distributions
dist_arr        = reshape(dist, glob.Nkf, glob.Nbf, glob.Nzf);
dist_pol_arr    = reshape(dist_pol, glob.Nkf, glob.Nbf, glob.Nzf);
figure;
subplot(1, 2, 1);
title('No policy');
surf(glob.kgridf, glob.bgridf, dist_arr(:, :, 1)');
xlabel('$K$');
ylabel('$D / K$')
zlabel('Density')
title(['Ergodic distribution, $Z = ', num2str(Zfix), '$'])

subplot(1, 2, 2);
title('Simple policy');
surf(glob.kgridf, glob.bgridf, dist_pol_arr(:, :, 1)');
xlabel('$K$');
ylabel('$D / K$')
zlabel('Density')
title(['Ergodic distribution, $Z = ', num2str(Zfix), '$'])

end