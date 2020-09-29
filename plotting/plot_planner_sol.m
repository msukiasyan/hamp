function plot_planner_sol(sol, eq, glob_eq, param, glob, options, Kfix, mufix, Zfix)
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


implied_tax_lab     = sol.implied_tax_lab;
implied_tax_r       = sol.implied_tax_r;
transfer            = sol.transfer;

%% Prepare objects
ps              = gridmake([Kfix * 0.9; Kfix; Kfix * 1.1], glob.bgridf, mufix, [Zfix * 0.9; Zfix; Zfix * 1.1]);

Phip            = funbas(glob.fspace, ps, [0, 0, 0, 0]);
tax_labf        = Phip * (Phil \ (Phiu \ implied_tax_lab));
tax_rf          = Phip * (Phil \ (Phiu \ implied_tax_r));
transferf       = Phip * (Phil \ (Phiu \ transfer));
Yf              = Phip * (Phil \ (Phiu \ sol.Y));

tax_labf_arr    = reshape(tax_labf, 3, glob.Nbf, 1, 3);
tax_rf_arr      = reshape(tax_rf, 3, glob.Nbf, 1, 3);
transferf_arr   = reshape(transferf, 3, glob.Nbf, 1, 3);
Yf_arr          = reshape(Yf, 3, glob.Nbf, 1, 3);

%% Plot taxes and transfers

% labor tax
figure;
subplot(1, 2, 1);
plot(glob.bgridf, tax_labf_arr(:, :, :, 1), 'LineWidth', 2);
yline(0, 'LineStyle', '--');
xlabel('$D / K$');
ylabel('$\tau_l$');
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Implied tax on labor, $Z = ' num2str(Zfix * 0.9) '$, $\mu = ', num2str(mufix), '$']);
grid on;

subplot(1, 2, 2);
plot(glob.bgridf, tax_labf_arr(:, :, :, 1) - tax_labf_arr(:, :, :, 2), 'LineWidth', 2);
yline(0, 'LineStyle', '--');
xlabel('$D / K$');
ylabel('$\Delta \tau_l$');
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Difference in $\tau_l$ between $Z = ' num2str(Zfix * 0.9) '$ and $Z$ at mean, $\mu = ', num2str(mufix), '$']);
grid on;

% r tax
figure;
subplot(1, 2, 1);
plot(glob.bgridf, tax_rf_arr(:, :, :, 1), 'LineWidth', 2);
yline(0, 'LineStyle', '--');
xlabel('$D / K$');
ylabel('$\tau_r$');
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Implied tax on dep. issuance, $Z = ' num2str(Zfix * 0.9) '$, $\mu = ', num2str(mufix), '$']);
grid on;

subplot(1, 2, 2);
plot(glob.bgridf, tax_rf_arr(:, :, :, 1) - tax_rf_arr(:, :, :, 2), 'LineWidth', 2);
yline(0, 'LineStyle', '--');
xlabel('$D / K$');
ylabel('$\Delta \tau_r$');
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Difference in $\tau_r$ between $Z = ' num2str(Zfix * 0.9) '$ and $Z$ at mean, $\mu = ', num2str(mufix), '$']);
grid on;

figure;
plot(glob.bgridf, transferf_arr(:, :, :, 1) ./ Yf_arr(:, :, :, 1), 'LineWidth', 2);
yline(0, 'LineStyle', '--');
xlabel('$D / K$');
ylabel('$T / Y$');
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Implied transfer as fraction of output, $Z = ' num2str(Zfix * 0.9) '$, $\mu = ', num2str(mufix), '$']);
grid on;

%% Plot equilibrium objects
ps          = gridmake([Kfix * 0.9; Kfix; Kfix * 1.1], glob.bgridf, 0, [Zfix * 0.9; Zfix; Zfix * 1.1]);

Phip        = funbas(glob.fspace, ps, [0, 0, 0, 0]);
Phip_eq     = funbas(glob_eq.fspace, ps(:, [1, 2, 4]), [0, 0, 0]);

I_pl        = Phip * (Phil \ (Phiu \ I));
I_eq        = Phip_eq * (glob_eq.Phil \ (glob_eq.Phiu \ eq.I));

I_pl_arr    = reshape(I_pl, 3, glob.Nbf, 3);
I_eq_arr    = reshape(I_eq, 3, glob.Nbf, 3);

L_pl        = Phip * (Phil \ (Phiu \ L));
L_eq        = Phip_eq * (glob_eq.Phil \ (glob_eq.Phiu \ eq.L));

L_pl_arr    = reshape(L_pl, 3, glob.Nbf, 3);
L_eq_arr    = reshape(L_eq, 3, glob.Nbf, 3);

c_w_pl      = Phip * (Phil \ (Phiu \ c_w));
c_w_eq      = Phip_eq * (glob_eq.Phil \ (glob_eq.Phiu \ eq.c_w));

c_w_pl_arr  = reshape(c_w_pl, 3, glob.Nbf, 3);
c_w_eq_arr  = reshape(c_w_eq, 3, glob.Nbf, 3);

c_b_pl      = Phip * (Phil \ (Phiu \ c_b));
c_b_eq      = Phip_eq * (glob_eq.Phil \ (glob_eq.Phiu \ eq.c_b));

c_b_pl_arr  = reshape(c_b_pl, 3, glob.Nbf, 3);
c_b_eq_arr  = reshape(c_b_eq, 3, glob.Nbf, 3);

figure;
subplot(2, 2, 1);
plot(glob.bgridf, c_w_pl_arr(:, :, 1), 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, c_w_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('$D / K$');
ylabel('$C_w$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Worker cons, $Z = ' num2str(Zfix * 0.9), '$'])
grid on

subplot(2, 2, 2);
plot(glob.bgridf, c_b_pl_arr(:, :, 1), 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, c_b_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('$D / K$');
ylabel('$C_e$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Expert cons, $Z = ' num2str(Zfix * 0.9), '$'])
grid on

subplot(2, 2, 3);
plot(glob.bgridf, I_pl_arr(:, :, 1), 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, I_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('$D / K$');
ylabel('$I$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Investment, $Z = ' num2str(Zfix * 0.9), '$'])
grid on

subplot(2, 2, 4);
plot(glob.bgridf, L_pl_arr(:, :, 1), 'LineWidth', 2)
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, L_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('$D / K$');
ylabel('$L$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Labor, $Z = ' num2str(Zfix * 0.9), '$'])
grid on

%% Plot welfare
ws                  = Phip * sol.cold;
omeg                = utility_to_pce(ws, param, glob, options);
omeg_arr            = reshape(omeg, 3, glob.Nbf, 3);

[c_v_w, c_v_b]      = compute_value_function(eq, param, glob_eq, options);
ws_eq               = Phip_eq * (glob.lambda * c_v_b + (1 - glob.lambda) * c_v_w);
omeg_eq             = utility_to_pce(ws_eq, param, glob, options);
omeg_eq_arr         = reshape(omeg_eq, 3, glob.Nbf, 3);

figure;
plot(glob.bgridf, omeg_arr(:, :, 1), 'LineWidth', 2);
hold on
set(gca,'ColorOrderIndex',1)
plot(glob.bgridf, omeg_eq_arr(:, :, 1), 'LineStyle', '--', 'LineWidth', 2)
xlabel('$D / K$');
ylabel('$\omega$');
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Planner''s objective, $Z = ' num2str(Zfix * 0.9) '$, $\mu = ', num2str(mufix), '$']);
grid on;

end