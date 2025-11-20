clear; clc;


output_dir = 'PCR_Results_Final';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

analyses = struct();

% Adversity Analysis Configuration
analyses(1).name = 'Adversity_Amygdala';
analyses(1).X = table2array(data_adv(:, 6:30));
analyses(1).X(analyses(1).X == 2) = 0; 
analyses(1).X(isnan(analyses(1).X)) = 0;
analyses(1).y = data_adv.deltar_amyg_zscore;
analyses(1).labels = adv_count.Properties.VariableNames(2:26);

% Cognitive Analysis Configuration
analyses(2).name = 'Cognition_Hippocampus';
analyses(2).X = table2array(cogcomp_table(:, 2:8));
analyses(2).y = cogcomp_table.deltar_h_zscore;
analyses(2).labels = {'tpvt','picseq','patterncomp','read','flanker','dccs','lswmt'};

n_folds = 10;
n_bootstraps = 10000; 

%% Main Analysis Loop
for a = 1:length(analyses)
    fprintf('Running Analysis: %s...\n', analyses(a).name);
    
    X = analyses(a).X;
    y = analyses(a).y;
    feat_labels = analyses(a).labels;
    [n_samples, n_features] = size(X);
    
    %% Cross-Validation for Optimal Components
    max_components = min(n_features, n_samples - 1); 
    mse_pcr = zeros(max_components, 1);
    cv = cvpartition(n_samples, 'KFold', n_folds);
    
    for n_comp = 1:max_components
        y_pred_fold = zeros(n_samples, 1);
        
        for i = 1:cv.NumTestSets
            tr_idx = cv.training(i);
            te_idx = cv.test(i);
            
            [X_tr, mu, sigma] = zscore(X(tr_idx, :));
            X_te = (X(te_idx, :) - mu) ./ sigma;
            X_te(:, sigma == 0) = 0; 
            
            [coeff_tr, score_tr] = pca(X_tr, 'NumComponents', n_comp);
            
            X_reg_tr = [ones(sum(tr_idx), 1), score_tr];
            b = X_reg_tr \ y(tr_idx);
            
            score_te = X_te * coeff_tr;
            X_reg_te = [ones(sum(te_idx), 1), score_te];
            y_pred_fold(te_idx) = X_reg_te * b;
        end
        mse_pcr(n_comp) = mean((y_pred_fold - y).^2);
    end
    
    [~, best_n_comp] = min(mse_pcr);
    
    %% Bootstrapping for Feature Importance
    boot_weights = zeros(n_features, n_bootstraps);
    
    for b = 1:n_bootstraps
        boot_idx = randi(n_samples, n_samples, 1);
        X_boot = X(boot_idx, :);
        y_boot = y(boot_idx);
        
        [X_boot_z, ~, sigma_b] = zscore(X_boot);
        X_boot_z(:, sigma_b == 0) = 0;
        
        [coeff_b, score_b] = pca(X_boot_z, 'NumComponents', best_n_comp);
        
        X_reg_b = [ones(n_samples, 1), score_b];
        beta_b = X_reg_b \ y_boot;
        
        boot_weights(:, b) = coeff_b * beta_b(2:end);
    end
    
    %% Statistics and File Saving
    mean_weights = mean(boot_weights, 2);
    ci_95 = prctile(boot_weights, [2.5, 97.5], 2);
    
    p_vals = zeros(n_features, 1);
    for f = 1:n_features
        if mean_weights(f) > 0
            p_vals(f) = 2 * mean(boot_weights(f, :) <= 0);
        else
            p_vals(f) = 2 * mean(boot_weights(f, :) >= 0);
        end
    end
    
    T = table(feat_labels', mean_weights, ci_95(:,1), ci_95(:,2), p_vals, ...
        'VariableNames', {'Feature', 'Importance', 'CI_Lower', 'CI_Upper', 'P_Value'});
    
    [~, sort_idx] = sort(abs(T.Importance), 'descend');
    T = T(sort_idx, :);
    T.Significant = (T.CI_Lower .* T.CI_Upper) > 0;
    
    writetable(T, fullfile(output_dir, [analyses(a).name '_Results.csv']));
    
    %% Plot Generation
    f = figure('Visible', 'off', 'Position', [100, 100, 1200, 600]);
    hold on;
    
    bar_data = T.Importance;
    err_low = T.Importance - T.CI_Lower;
    err_high = T.CI_Upper - T.Importance;
    
    b_plot = bar(bar_data);
    b_plot.FaceColor = [0.4 0.6 0.8];
    
    errorbar(1:n_features, bar_data, err_low, err_high, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    
    title(strrep(analyses(a).name, '_', ' '), 'FontSize', 14);
    ylabel('Predictive Importance');
    xticks(1:n_features);
    xticklabels(T.Feature);
    xtickangle(45);
    grid on;
    
    sig_idx = find(T.Significant);
    if ~isempty(sig_idx)
        y_limits = ylim;
        offset = (y_limits(2) - y_limits(1)) * 0.05;
        plot(sig_idx, T.CI_Upper(sig_idx) + offset, 'r*', 'MarkerSize', 8);
    end
    
    saveas(f, fullfile(output_dir, [analyses(a).name '_FeatureImportance.png']));
    close(f);
end

disp('Analysis complete.');