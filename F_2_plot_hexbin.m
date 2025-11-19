%% Adversity
x = 50 - hcpd_data.ale_total_number;
y = hcpd_data.deltar_a_zscore;

figure('Color', 'w');
ax = gca;
hold(ax, 'on');

hexscatter(x, y, 'res', 27);
colormap(ax, 'parula'); 

mdl = fitlm(x, y);
x_range = linspace(min(x), max(x), 100)';
[y_pred, y_ci] = predict(mdl, x_range);

fill(ax, [x_range; flipud(x_range)], [y_ci(:,1); flipud(y_ci(:,2))], 'r', ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(ax, x_range, y_pred, 'r-', 'LineWidth', 2);

hold(ax, 'off');
box off;
set(gca, 'TickDir', 'out');
xlabel('Cumulative Adversity');
ylabel('Amygdala delta-r Z-score');
title('Adversity vs Amygdala Differentiation');

%% Cognition
x = cog_deltar_z.nih_totalcogcomp_ageadjusted;
y = cog_deltar_z.deltar_amyg_zscore; 

figure('Color', 'w');
ax = gca;
hold(ax, 'on');

hexscatter(x, y, 'res', 23);
colormap(ax, 'parula'); 

mdl = fitlm(x, y);
x_range = linspace(min(x), max(x), 100)';
[y_pred, y_ci] = predict(mdl, x_range);

fill(ax, [x_range; flipud(x_range)], [y_ci(:,1); flipud(y_ci(:,2))], 'r', ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(ax, x_range, y_pred, 'r-', 'LineWidth', 2);

hold(ax, 'off');
box off;
set(gca, 'TickDir', 'out', 'TickLength', [0 0]);
xlabel('Cognitive Ability (Age Adjusted)');
ylabel('Amygdala delta-r Z-score'); 
title('Cognition vs Amygdala Differentiation');