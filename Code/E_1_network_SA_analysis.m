clear; clc;
addpath('/Users/wil/research/tools/');

dir_amyg = 'SAaxis_network_mapping_amyg_detailed';
dir_hipp = 'SAaxis_network_mapping_hipp_detailed';

file_glasser = 'Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii';
file_schaefer = 'Schaefer2018_400Parcels_7Networks_order_Tian_Subcortex_S1.dlabel.nii';
file_sa_axis = 'SAaxis_cifti_label_frontal.csv';

ages = (0.5:1:21.5)';
frontal_ids = [42,57,59,73,74,75,76,77,78,79,80,81,85,86,88,104,106,109,180,181,182,195,196,222,237,239,253,254,255,256,257,258,259,260,261,265,266,268,284,286,289,360,361,362,375,376,26,27,28,60,84,87,93,94,96,97,98,100,101,102,103,105,107,112,113,186,187,206,207,208,240,264,267,273,274,277,278,280,281,282,283,285,287,292,293,366,367,82,83,89,90,91,92,95,99,108,110,114,262,263,269,270,271,272,275,276,279,288,290,294];
net_ids = [2, 3, 4, 5, 6, 7];
net_names = {'SomMot','DorsAttn','Salience','Control','Default','Limbic'};
net_colors = [70 130 180; 0 117 16; 196 58 251; 232 147 33; 205 61 80; 219 248 165] / 255;

%% Load Atlases and Create Masks
glasser = cifti_read(file_glasser);
schaefer = cifti_read(file_schaefer);
sa_raw = readmatrix(file_sa_axis);

is_frontal = false(size(glasser.cdata));
is_frontal(ismember(glasser.cdata, frontal_ids)) = true;

labels = schaefer.diminfo{1,2}.maps.table; 
map_net = nan(size(schaefer.cdata));
map_hemi = nan(size(schaefer.cdata));

for i = 1:length(labels)
    name = lower(labels(i).name);
    key = labels(i).key;
    idx = (schaefer.cdata == key);
    
    id = 8; 
    if contains(name, 'vis'), id=1; elseif contains(name, 'sommot'), id=2;
    elseif contains(name, 'dorsattn'), id=3; elseif contains(name, 'salience'), id=4;
    elseif contains(name, 'control'), id=5; elseif contains(name, 'default'), id=6;
    elseif contains(name, 'limbic'), id=7; end
    map_net(idx) = id;
    
    if i >= 17 && i <= 216, map_hemi(idx) = 1;
    elseif i >= 217 && i <= 416, map_hemi(idx) = 2;
    end
end

map_sa = nan(size(glasser.cdata));
if length(sa_raw) == sum(is_frontal)
    map_sa(is_frontal) = sa_raw;
else
    map_sa = sa_raw;
end

mask_valid = is_frontal & ~isnan(map_net) & ~isnan(map_hemi) & ~isnan(map_sa);
data_net = map_net(mask_valid);
data_hemi = map_hemi(mask_valid);
data_sa = map_sa(mask_valid);

%% Load Connectivity Data and Organize
Y_vec = []; X_age = []; X_sa = []; X_circ = []; X_net = [];

nBands = length(ages);
nGroups = length(net_ids) * 2; 
Store_Means_H = nan(nBands, nGroups);
Store_Means_A = nan(nBands, nGroups);
Store_SA_Pos  = nan(nGroups, 1);
Group_Labels  = cell(nGroups, 1);

dirs = {dir_hipp, dir_amyg};

g_count = 0;
for n = 1:length(net_ids)
    for h = 1:2
        g_count = g_count + 1;
        mask_nh = (data_net == net_ids(n)) & (data_hemi == h);
        Store_SA_Pos(g_count) = mean(data_sa(mask_nh), 'omitnan');
        hemi_str = 'LH'; if h==2, hemi_str='RH'; end
        Group_Labels{g_count} = [hemi_str '-' net_names{n}];
    end
end

for c = 1:2
    files = dir(fullfile(dirs{c}, '*.dscalar.nii'));
    [~, order] = sort({files.name});
    files = files(order);
    
    for b = 1:length(files)
        cii = cifti_read(fullfile(files(b).folder, files(b).name));
        vals = cii.cdata(mask_valid);
        vals_z = (vals - mean(vals)) / std(vals);
        
        g_count = 0;
        for n = 1:length(net_ids)
            nid = net_ids(n);
            for h = 1:2
                g_count = g_count + 1;
                mask_nh = (data_net == nid) & (data_hemi == h);
                
                avg_conn = mean(vals_z(mask_nh), 'omitnan');
                sa_val   = mean(data_sa(mask_nh), 'omitnan');
                
                Y_vec = [Y_vec; avg_conn];
                X_age = [X_age; ages(b)];
                X_sa  = [X_sa; sa_val];
                X_circ= [X_circ; (c-1)]; 
                X_net = [X_net; nid];
                
                if c == 1
                    Store_Means_H(b, g_count) = avg_conn;
                else
                    Store_Means_A(b, g_count) = avg_conn;
                end
            end
        end
    end
end

%% Stats
Z_age = (X_age - mean(X_age)) / std(X_age);
Z_sa  = (X_sa - mean(X_sa)) / std(X_sa);

X_mat = [ones(size(Y_vec)), Z_age, Z_sa, Z_age.*Z_sa, X_circ, ...
         Z_age.*X_circ, Z_sa.*X_circ, Z_age.*Z_sa.*X_circ];

[b_glob, ~, ~, ~, ~] = regress(Y_vec, X_mat);

mask_sens = (X_net ~= 4) & (X_net ~= 6);
[b_sens, ~, ~, ~, ~] = regress(Y_vec(mask_sens), X_mat(mask_sens,:));

Diff_Traj = Store_Means_A - Store_Means_H;
Slopes = nan(nGroups, 1);
X_design = [ones(nBands, 1), ages];

for g = 1:nGroups
    if all(~isnan(Diff_Traj(:,g)))
        beta = X_design \ Diff_Traj(:,g);
        Slopes(g) = beta(2);
    end
end

[rho, pval] = corr(Store_SA_Pos, Slopes, 'Type', 'Spearman');

%% Visualization
figure('Color','w', 'Position', [100 100 1200 500]); tiledlayout(1, 3);

titles = {'Hippocampus Model', 'Amygdala Model'};
x_sim = linspace(min(Z_age), max(Z_age), 50)';
sa_levels = [-1, 1];
colors_sa = [0.2 0.2 0.2; 0.6 0.6 0.6];

nexttile; hold on; title('Interaction: Age \times SA \times Circuit');
calc_pred = @(age, sa, c) b_glob(1) + b_glob(2)*age + b_glob(3)*sa + ...
    b_glob(4)*(age.*sa) + b_glob(5)*c + b_glob(6)*(age.*c) + ...
    b_glob(7)*(sa.*c) + b_glob(8)*(age.*sa.*c);

for i = 1:2
    plot(x_sim, calc_pred(x_sim, sa_levels(i), 0), '-', 'LineWidth', 2, 'Color', colors_sa(i,:));
end
for i = 1:2
    plot(x_sim, calc_pred(x_sim, sa_levels(i), 1), '--', 'LineWidth', 2, 'Color', colors_sa(i,:));
end
legend('Hipp (Sens)','Hipp (Assoc)','Amyg (Sens)','Amyg (Assoc)','Location','best');

nexttile; hold on;
scatter(Store_SA_Pos, Slopes, 60, 'filled', 'k');
lsline;
for g = 1:nGroups
    text(Store_SA_Pos(g)+100, Slopes(g), Group_Labels{g}, 'FontSize', 7, 'Interpreter','none');
end
title(sprintf('Differentiation vs Hierarchy\n\\rho = %.3f', rho));

nexttile; hold on; title('Hippocampus Network Trajectories');
for i = 1:nBands
    if mod(i,2)==0 
        patch([ages(i)-0.5 ages(i)+0.5 ages(i)+0.5 ages(i)-0.5], ...
              [-1 -1 1 1], [0.9 0.9 0.9], 'EdgeColor', 'none');
    end
end

for n = 1:length(net_ids)
    cols_n = contains(Group_Labels, net_names{n});
    traj = mean(Store_Means_H(:, cols_n), 2);
    plot(ages, traj, '-o', 'Color', net_colors(n,:), ...
        'LineWidth', 2, 'MarkerFaceColor', net_colors(n,:));
end
ylim([-0.6 0.6]);