function plot_eqbm(eq, param, glob, options, Kfix, Zfix)
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
s           = glob.s;
Phi         = glob.Phi;
Phiu        = glob.Phiu;
Phil        = glob.Phil;

%% Plot
ps          = gridmake([Kfix * 0.9; Kfix; Kfix * 1.1], glob.bgridf, Zfix);

Phip        = funbas(glob.fspace, ps);

If          = Phip * (Phil \ (Phiu \ I));
If_arr      = reshape(If, 3, glob.Nbf);

Lf          = Phip * (Phil \ (Phiu \ L));
Lf_arr      = reshape(Lf, 3, glob.Nbf);

rf          = Phip * (Phil \ (Phiu \ r));
rf_arr      = reshape(rf, 3, glob.Nbf);

c_wf        = Phip * (Phil \ (Phiu \ c_w));
c_wf_arr    = reshape(c_wf, 3, glob.Nbf);

c_bf        = Phip * (Phil \ (Phiu \ c_b));
c_bf_arr    = reshape(c_bf, 3, glob.Nbf);

figure;

subplot(1, 3, 1);
plot(glob.bgridf, If_arr, 'LineWidth', 2)
xlabel('$D / K$')
ylabel('$I$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Investment, $Z = ', num2str(Zfix), '$'])
grid on

subplot(1, 3, 2);
plot(glob.bgridf, Lf_arr, 'LineWidth', 2)
xlabel('$D / K$')
ylabel('$L$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Labor, $Z = ', num2str(Zfix), '$'])
grid on

subplot(1, 3, 3);
plot(glob.bgridf, rf_arr, 'LineWidth', 2)
xlabel('$D / K$')
ylabel('$r$')
legend(['$K$ at -10\%'], ['$K$ at mean'], ['$K$ at +10\%']);
title(['Deposit rate, $Z = ', num2str(Zfix), '$'])
grid on

dist_arr = reshape(dist, glob.Nkf, glob.Nbf, glob.Nzf);
figure;
surf(glob.kgridf, glob.bgridf, dist_arr(:, :, 1)');
xlabel('$K$');
ylabel('$D / K$')
zlabel('Density')
title('Ergodic distribution, $Z$ at low')
end