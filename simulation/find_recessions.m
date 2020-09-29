function reces = find_recessions(sim, param, glob, options, reces_times)
%% Unpack
burnin          = options.burnin + 1;
expi            = 1;
Kt              = sim.Kt(expi, burnin:end)';
Bt              = sim.Bt(expi, burnin:end)';
Zt              = sim.Zt(expi, burnin:end)';
c_wt            = sim.c_wt(expi, burnin:end)';
c_bt            = sim.c_bt(expi, burnin:end)';
Yt              = sim.Yt(expi, burnin:end)';
Lt              = sim.Lt(expi, burnin:end)';
It              = sim.It(expi, burnin:end)';
qt              = sim.qt(expi, burnin:end)';
rt              = sim.rt(expi, burnin:end)';
levt            = sim.levt(expi, burnin:end)';
T               = size(Kt, 1);

%% Find recession episodes
if nargin < 5
    stdB            = std(Bt);
    meanB           = mean(Bt);
    ratt            = Bt ./ Kt;
    
    stdrat          = std(ratt);
    meanrat         = mean(ratt);

    stdlev          = std(levt);
    meanlev         = mean(levt);

    Nrec            = 0;
    reces_times     = zeros(T, 1);

    for t = 41:T-20
        if levt(t) > meanlev + 2 * stdlev && levt(t - 1) <= meanlev + 2 * stdlev && Zt(t) == glob.zgrid(1)
        % if ratt(t) > meanrat + 2 * stdrat && ratt(t - 1) <= meanrat + 2 * stdrat && Zt(t) == glob.zgrid(1)
            reces_times(Nrec + 1)   = t;
            Nrec    = Nrec + 1;
        end
    end
    reces_times     = reces_times(1:Nrec);
end

Nrec            = size(reces_times, 1);

recKt           = zeros(Nrec, 61);
recBt           = zeros(Nrec, 61);
recZt           = zeros(Nrec, 61);
recc_wt         = zeros(Nrec, 61);
recc_bt         = zeros(Nrec, 61);
recYt           = zeros(Nrec, 61);
recLt           = zeros(Nrec, 61);
recIt           = zeros(Nrec, 61);
recqt           = zeros(Nrec, 61);
recrt           = zeros(Nrec, 61);
reclevt         = zeros(Nrec, 61);

for ri = 1:Nrec
    recKt(ri, :)        = Kt(reces_times(ri) - 40:reces_times(ri) + 20);
    recBt(ri, :)        = Bt(reces_times(ri) - 40:reces_times(ri) + 20);
    recZt(ri, :)        = Zt(reces_times(ri) - 40:reces_times(ri) + 20);
    recc_wt(ri, :)      = c_wt(reces_times(ri) - 40:reces_times(ri) + 20);
    recc_bt(ri, :)      = c_bt(reces_times(ri) - 40:reces_times(ri) + 20);
    recYt(ri, :)        = Yt(reces_times(ri) - 40:reces_times(ri) + 20);
    recLt(ri, :)        = Lt(reces_times(ri) - 40:reces_times(ri) + 20);
    recIt(ri, :)        = It(reces_times(ri) - 40:reces_times(ri) + 20);
    recqt(ri, :)        = qt(reces_times(ri) - 40:reces_times(ri) + 20);
    recrt(ri, :)        = rt(reces_times(ri) - 40:reces_times(ri) + 20);
    reclevt(ri, :)      = levt(reces_times(ri) - 40:reces_times(ri) + 20);
end

Ktmean          = mean(recKt, 1);
Btmean          = mean(recBt, 1);
Ztmean          = mean(recZt, 1);
c_wtmean        = mean(recc_wt, 1);
c_btmean        = mean(recc_bt, 1);
Ytmean          = mean(recYt, 1);
Ltmean          = mean(recLt, 1);
Itmean          = mean(recIt, 1);
qtmean          = mean(recqt, 1);
rtmean          = mean(recrt, 1);
levtmean        = mean(reclevt, 1);

%% Pack
reces.reces_times   = reces_times;
reces.Ktmean        = Ktmean;
reces.Btmean        = Btmean;
reces.Ztmean        = Ztmean;
reces.c_wtmean      = c_wtmean;
reces.c_btmean      = c_btmean;
reces.Ytmean        = Ytmean;
reces.Ltmean        = Ltmean;
reces.Itmean        = Itmean;
reces.qtmean        = qtmean;
reces.rtmean        = rtmean;
reces.levtmean      = levtmean;

end