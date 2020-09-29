function plot_welfare_simple(ceq_w, ceq_b, cap, shock, param, glob, options)
%% Plot
ceq_w_arr       = reshape(ceq_w * 100 - 100, glob.Nk, glob.Nb, glob.Nz);
ceq_b_arr       = reshape(ceq_b * 100 - 100, glob.Nk, glob.Nb, glob.Nz);

finegrid        = linspace(min(glob.bgrid), max(glob.bgrid), 100);

ceq_w_coef      = funfitxy(glob.fspace, glob.s, ceq_w);
ceq_w_val       = funeval(ceq_w_coef, glob.fspace, [repmat(cap, 100, 1), finegrid', repmat(shock, 100, 1)]);

ceq_b_coef      = funfitxy(glob.fspace, glob.s, ceq_b);
ceq_b_val       = funeval(ceq_b_coef, glob.fspace, [repmat(cap, 100, 1), finegrid', repmat(shock, 100, 1)]);

figure;
subplot(1, 2, 1)
plot(finegrid, ceq_w_val * 100 - 100, 'LineWidth', 2);
xline(glob.target_lev, '--k');
xlabel('D / K')
ylabel('% Consumption gain')
title('Welfare gain from simple policy, workers')
grid on

subplot(1, 2, 2)
plot(finegrid, ceq_b_val * 100 - 100, 'LineWidth', 2);
xline(glob.target_lev, '--k');
xlabel('D / K')
ylabel('% Consumption gain')
title('Welfare gain from simple policy, experts')
grid on

end