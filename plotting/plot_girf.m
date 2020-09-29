function plot_girf(g_eq, g_plan, param, glob, options)
%% Plot

figure;
% sgtitle('Generalized IRFs at the ergodic mean of comp. eqbm.');

subplot(3, 2, 1);
plot(1:options.Tgirf, mean(g_eq.Zt, 1), 'LineWidth', 2)
yline(glob.zgrid' * glob.Pssz, 'LineStyle', '--');
xlabel('Time')
ylabel('TFP')
title('TFP')
grid on

subplot(3, 2, 2);
plot(1:options.Tgirf, [mean(g_eq.Yt, 1); mean(g_plan.Yt, 1)] ./ mean(g_eq.Yt(:, end)), 'LineWidth', 2)
yline(1, 'LineStyle', '--');
xlabel('Time')
ylabel('Output')
title('Output rel. to ergodic mean')
grid on

subplot(3, 2, 3);
plot(1:options.Tgirf, [mean(g_eq.It, 1); mean(g_plan.It, 1)] ./ mean(g_eq.It(:, end)), 'LineWidth', 2)
yline(1, 'LineStyle', '--');
xlabel('Time')
ylabel('Investment')
title('Investment rel. to ergodic mean')
grid on

subplot(3, 2, 4);
plot(1:options.Tgirf, [mean(g_eq.Lt, 1); mean(g_plan.Lt, 1)] ./ mean(g_eq.Lt(:, end)), 'LineWidth', 2)
yline(1, 'LineStyle', '--');
xlabel('Time')
ylabel('Labor')
title('Labor rel. to ergodic mean')
grid on

subplot(3, 2, 5);
plot(1:options.Tgirf, [mean(g_eq.Kt, 1); mean(g_plan.Kt, 1)] ./ mean(g_eq.Kt(:, end)), 'LineWidth', 2)
yline(1, 'LineStyle', '--');
xlabel('Time')
ylabel('Capital')
title('Capital rel. to ergodic mean')
grid on

subplot(3, 2, 6);
plot(1:options.Tgirf, [mean(g_eq.rt, 1); mean(g_plan.rt, 1)], 'LineWidth', 2)
yline(g_eq.r0, 'LineStyle', '--');
xlabel('Time')
ylabel('Risk-free rate')
title('Risk-free rate')
grid on

figure;
sgtitle('Generalized IRFs, taxes/transfers');

subplot(1, 3, 1);
plot(1:options.Tgirf, mean(g_plan.tax_labt, 1) * 100, 'LineWidth', 2)
yline(0, 'LineStyle', '--');
xlabel('Time')
ylabel('Tax, \%')
title('Labor tax, \%')
grid on

subplot(1, 3, 2);
plot(1:options.Tgirf, mean(g_plan.tax_rt, 1) * 100, 'LineWidth', 2)
yline(0, 'LineStyle', '--');
xlabel('Time')
ylabel('Tax, \%')
title('Deposit iss. tax, \%')
grid on

subplot(1, 3, 3);
plot(1:options.Tgirf, mean(g_plan.transft ./ g_plan.Yt, 1) * 100, 'LineWidth', 2)
yline(0, 'LineStyle', '--');
xlabel('Time')
ylabel('$T/Y$, \%')
title('Transfer as \% of output')
grid on

% subplot(3, 2, 5);
% plot(1:options.Tgirf, c_wtmean ./ mean(sim.c_wt, 2), 'LineWidth', 2)
% xlabel('Time')
% ylabel('Worker consumption')
% title('Worker consumption rel. to ergodic mean')
% grid on
% 
% subplot(3, 2, 6);
% plot(1:options.Tgirf, c_btmean ./ mean(sim.c_bt, 2), 'LineWidth', 2)
% xlabel('Time')
% ylabel('Expert consumption')
% title('Expert consumption rel. to ergodic mean')
% grid on

end