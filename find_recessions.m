function reces = find_recessions(sim, param, glob, options)
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
stdB            = std(Bt);
meanB           = mean(Bt);

Nrec            = 0;
reces_times     = zeros(T, 1);

for t = 41:T-20
    if Bt(t) < meanB - 2 * stdB && Bt(t - 1) >= meanB - 2 * stdB
        reces_times(Nrec + 1)   = t;
        Nrec    = Nrec + 1;
    end
end
reces_times     = reces_times(1:Nrec);
% reces_times     = find(Bt < meanB - 2 * stdB);
% Nrec            = size(reces_times, 1);
% 
% recNu           = Nrec;
% recNl           = 1;
% 
% for ri = Nrec:(-1):1
%     if reces_times(ri) + 20 > T
%         recNu = ri - 1;
%     end
%     if reces_times(ri) - 40 < 1
%         if reces_times(ri + 1) - 40 >= 1
%             recNl = ri + 1;
%         end
%     end
% end
% reces_times     = reces_times(recNl:recNu);
% Nrec            = size(reces_times, 1);

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

end