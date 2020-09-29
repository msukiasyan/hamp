function plot_recessions(reces, sim, param, glob, options, reces_pol)
%% Unpack
reces_times     = reces.reces_times;
Ktmean          = reces.Ktmean;
Btmean          = reces.Btmean;
Ztmean          = reces.Ztmean;
c_wtmean        = reces.c_wtmean;
c_btmean        = reces.c_btmean;
Ytmean          = reces.Ytmean;
Ltmean          = reces.Ltmean;
Itmean          = reces.Itmean;
qtmean          = reces.qtmean;
rtmean          = reces.rtmean;
levtmean        = reces.levtmean;

if nargin > 5
    Ktmean          = [Ktmean; reces_pol.Ktmean];
    Btmean          = [Btmean; reces_pol.Btmean];
    Ztmean          = [Ztmean; reces_pol.Ztmean];
    c_wtmean        = [c_wtmean; reces_pol.c_wtmean];
    c_btmean        = [c_btmean; reces_pol.c_btmean];
    Ytmean          = [Ytmean; reces_pol.Ytmean];
    Ltmean          = [Ltmean; reces_pol.Ltmean];
    Itmean          = [Itmean; reces_pol.Itmean];
    qtmean          = [qtmean; reces_pol.qtmean];
    rtmean          = [rtmean; reces_pol.rtmean];
    levtmean        = [levtmean; reces_pol.levtmean];
end

%% Plot

figure;
sgtitle('Variables over crises');

subplot(3, 2, 1);
plot(-40:20, Itmean ./ mean(sim.It, 2), 'LineWidth', 2)
xlabel('Time')
ylabel('Mean Investment')
title('Investment rel. to ergodic mean')
grid on

subplot(3, 2, 2);
plot(-40:20, Ytmean ./ mean(sim.Yt), 'LineWidth', 2)
xlabel('Time')
ylabel('Output')
title('Output rel. to ergodic mean')
grid on

subplot(3, 2, 3);
plot(-40:20, rtmean, 'LineWidth', 2)
xlabel('Time')
ylabel('Risk-free rate')
title('Risk-free rate')
grid on

subplot(3, 2, 4);
plot(-40:20, Ztmean, 'LineWidth', 2)
xlabel('Time')
ylabel('TFP')
title('TFP')
grid on

subplot(3, 2, 5);
plot(-40:20, c_wtmean ./ mean(sim.c_wt, 2), 'LineWidth', 2)
xlabel('Time')
ylabel('Worker consumption')
title('Worker consumption rel. to ergodic mean')
grid on

subplot(3, 2, 6);
plot(-40:20, c_btmean ./ mean(sim.c_bt, 2), 'LineWidth', 2)
xlabel('Time')
ylabel('Expert consumption')
title('Expert consumption rel. to ergodic mean')
grid on

end