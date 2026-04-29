clear; clc;
Top = '/data/projects/punim2400';
addpath(genpath(fullfile(Top, 'tools')));
out_dir = fullfile(Top, 'derivatives', 'voxelwise_DI_age_effect');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
cov_table = readtable(fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_table_BCPHCPD.csv'));
delta_r_voxels = readmatrix(fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_voxel_BCPHCPD.csv'));
age = cov_table.interview_age;
n_vox_common = size(delta_r_voxels, 2);
beta_vals = zeros(n_vox_common, 1);
for v = 1:n_vox_common
    mdl = fitlm(age, delta_r_voxels(:, v));
    beta_vals(v) = mdl.Coefficients.Estimate(2);
end
parc = cifti_label.cdata;
subcort_all_idx = find(ismember(parc, [1, 2, 9, 10]));
betamap = cifti_label;
betamap.cdata = zeros(size(parc, 1), 1);
betamap.cdata(subcort_all_idx) = beta_vals;
betamap.diminfo{2} = cifti_diminfo_make_scalars(1);
cifti_write(betamap, fullfile(out_dir, 'Age_Effect_Betamap.dscalar.nii'));
%% Stats & Predictions
n_vox = length(subcort_all_idx);
models = cifti_label.diminfo{1}.models(3:21);
ijk_map = containers.Map('KeyType','double','ValueType','any');
for m = 1:length(models)
    model = models{m};
    start_idx = model.start;
    count = model.count;
    vox_coords = model.voxlist;
    for i = 1:count
        ijk_map(start_idx + i - 1) = vox_coords(:, i);
    end
end
coords_mat = zeros(n_vox, 3);
for i = 1:n_vox
    ijk = ijk_map(subcort_all_idx(i));
    ijk_homogeneous = [ijk; 1];
    xyz_homogeneous = cifti_label.diminfo{1,1}.vol.sform * ijk_homogeneous;
    coords_mat(i, :) = xyz_homogeneous(1:3)';
end
subcort_betas = betamap.cdata(subcort_all_idx);
subcort_parc = parc(subcort_all_idx);
roi_idx = [1, 9, 2, 10];
all_preds_raw = cell(1, 4);
all_preds_norm = cell(1, 4);
all_xyz = cell(1, 4);
for r = 1:4
    idx_roi = subcort_parc == roi_idx(r);
    xyz_raw = coords_mat(idx_roi, :);
    betas = subcort_betas(idx_roi);
    
    xyz_model = xyz_raw;
    xyz_model(:, 1) = abs(xyz_model(:, 1));
    
    mdl = fitlm(xyz_model, betas, 'PredictorVars', {'X', 'Y', 'Z'});
    preds = predict(mdl, xyz_model);
    
    all_preds_raw{r} = preds;
    all_preds_norm{r} = 2 * ((preds - min(preds)) / (max(preds) - min(preds))) - 1;
    all_xyz{r} = xyz_raw;
end
% Global normalisation bounds
global_preds_raw_mat = cell2mat(all_preds_raw');
g_min = min(global_preds_raw_mat);
g_max = max(global_preds_raw_mat);
all_preds_global_norm = cell(1, 4);
for r = 1:4
    all_preds_global_norm{r} = 2 * ((all_preds_raw{r} - g_min) / (g_max - g_min)) - 1;
end
% Write Normalised CIFTIs concordant with plots
cdata_global = zeros(size(parc, 1), 1);
cdata_hippo = zeros(size(parc, 1), 1);
cdata_amyg = zeros(size(parc, 1), 1);
for r = 1:4
    idx_roi = subcort_parc == roi_idx(r);
    full_idx = subcort_all_idx(idx_roi);
    
    cdata_global(full_idx) = all_preds_global_norm{r};
    if r == 1 || r == 2
        cdata_hippo(full_idx) = all_preds_norm{r};
    elseif r == 3 || r == 4
        cdata_amyg(full_idx) = all_preds_norm{r};
    end
end
predmap = cifti_label;
predmap.diminfo{2} = cifti_diminfo_make_scalars(1);
predmap.cdata = cdata_global;
cifti_write(predmap, fullfile(out_dir, 'Gradient_Global.dscalar.nii'));
predmap.cdata = cdata_hippo;
cifti_write(predmap, fullfile(out_dir, 'Gradient_Hippocampus.dscalar.nii'));
predmap.cdata = cdata_amyg;
cifti_write(predmap, fullfile(out_dir, 'Gradient_Amygdala.dscalar.nii'));
% Run Bilateral Unnormalised Stats for Export
model_stats_table = table();
bilat_idx = {[1, 9], [2, 10]};
bilat_names = {'Hippocampus', 'Amygdala'};
for b = 1:2
    idx_roi = subcort_parc == bilat_idx{b}(1) | subcort_parc == bilat_idx{b}(2);
    xyz_raw = coords_mat(idx_roi, :);
    betas = subcort_betas(idx_roi);
    
    xyz_model = xyz_raw;
    xyz_model(:, 1) = abs(xyz_model(:, 1));
    
    mdl = fitlm(xyz_model, betas, 'PredictorVars', {'X', 'Y', 'Z'});
    
    coeff_tbl = mdl.Coefficients;
    Term = coeff_tbl.Properties.RowNames;
    ROI = repmat({bilat_names{b}}, length(Term), 1);
    temp_tbl = table(ROI, Term);
    temp_tbl = [temp_tbl, coeff_tbl];
    temp_tbl.Properties.RowNames = {};
    model_stats_table = [model_stats_table; temp_tbl];
end
writetable(model_stats_table, fullfile(out_dir, 'Spatial_Gradients_Model_Output.csv'));
%% Plotting
global_xyz = cell2mat(all_xyz'); 
x_b = [min(global_xyz(:,1)) - 5, max(global_xyz(:,1)) + 5];
y_b = [min(global_xyz(:,2)) - 5, max(global_xyz(:,2)) + 5];
z_b = [min(global_xyz(:,3)) - 5, max(global_xyz(:,3)) + 5];
[cx, cy, cz] = meshgrid(x_b, y_b, z_b);
corners_x = cx(:); corners_y = cy(:); corners_z = cz(:);
az = -115; 
el = 30;
c_anchors = [0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0];
custom_cmap = interp1(linspace(0,1,5), c_anchors, linspace(0,1,256));

% Global plot
fig1 = figure('Color', 'w', 'Position', [100, 100, 900, 700]);
hold on;
plot3([x_b(1) x_b(2) x_b(2) x_b(1) x_b(1)], [y_b(1) y_b(1) y_b(2) y_b(2) y_b(1)], [z_b(1) z_b(1) z_b(1) z_b(1) z_b(1)], 'k-', 'LineWidth', 0.5);
plot3([x_b(2) x_b(2) x_b(2) x_b(2) x_b(2)], [y_b(1) y_b(2) y_b(2) y_b(1) y_b(1)], [z_b(1) z_b(1) z_b(2) z_b(2) z_b(1)], 'k-', 'LineWidth', 0.5);
plot3([x_b(1) x_b(2) x_b(2) x_b(1) x_b(1)], [y_b(1) y_b(1) y_b(1) y_b(1) y_b(1)], [z_b(1) z_b(1) z_b(2) z_b(2) z_b(1)], 'k-', 'LineWidth', 0.5);

scatter3(corners_x, corners_y, corners_z, 1, 'w', 'MarkerEdgeAlpha', 0);
for r = 1:4
    xyz = all_xyz{r};
    preds = all_preds_global_norm{r};
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, preds, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
end
colormap(custom_cmap); 
clim([-1, 1]); 
c1 = colorbar;
axis equal; grid off;
xlim(x_b); ylim(y_b); zlim(z_b);
set(gca, 'PlotBoxAspectRatio', [1 1 1], 'DataAspectRatio', [1 1 1]);
set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
view([az, el]);
cam_pos = get(gca, 'CameraPosition');
cam_tgt = get(gca, 'CameraTarget');
cam_va  = get(gca, 'CameraViewAngle');
hold off;
saveas(fig1, fullfile(out_dir, 'Gradient_Global.png'));

% Amygdala plot
fig2 = figure('Color', 'w', 'Position', [150, 150, 900, 700]);
hold on;
plot3([x_b(1) x_b(2) x_b(2) x_b(1) x_b(1)], [y_b(1) y_b(1) y_b(2) y_b(2) y_b(1)], [z_b(1) z_b(1) z_b(1) z_b(1) z_b(1)], 'k-', 'LineWidth', 0.5);
plot3([x_b(2) x_b(2) x_b(2) x_b(2) x_b(2)], [y_b(1) y_b(2) y_b(2) y_b(1) y_b(1)], [z_b(1) z_b(1) z_b(2) z_b(2) z_b(1)], 'k-', 'LineWidth', 0.5);
plot3([x_b(1) x_b(2) x_b(2) x_b(1) x_b(1)], [y_b(1) y_b(1) y_b(1) y_b(1) y_b(1)], [z_b(1) z_b(1) z_b(2) z_b(2) z_b(1)], 'k-', 'LineWidth', 0.5);

scatter3(corners_x, corners_y, corners_z, 1, 'w', 'MarkerEdgeAlpha', 0);
for r = 3:4
    xyz = all_xyz{r};
    preds = all_preds_norm{r};
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, preds, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
end
colormap(custom_cmap); 
clim([-1, 1]); 
c2 = colorbar;
axis equal; grid off;
xlim(x_b); ylim(y_b); zlim(z_b);
set(gca, 'PlotBoxAspectRatio', [1 1 1], 'DataAspectRatio', [1 1 1]);
set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
set(gca, 'CameraPosition', cam_pos, 'CameraTarget', cam_tgt, 'CameraViewAngle', cam_va);
hold off;
saveas(fig2, fullfile(out_dir, 'Gradient_Amygdala.png'));

% Hippocampus plot
fig3 = figure('Color', 'w', 'Position', [200, 200, 900, 700]);
hold on;
plot3([x_b(1) x_b(2) x_b(2) x_b(1) x_b(1)], [y_b(1) y_b(1) y_b(2) y_b(2) y_b(1)], [z_b(1) z_b(1) z_b(1) z_b(1) z_b(1)], 'k-', 'LineWidth', 0.5);
plot3([x_b(2) x_b(2) x_b(2) x_b(2) x_b(2)], [y_b(1) y_b(2) y_b(2) y_b(1) y_b(1)], [z_b(1) z_b(1) z_b(2) z_b(2) z_b(1)], 'k-', 'LineWidth', 0.5);
plot3([x_b(1) x_b(2) x_b(2) x_b(1) x_b(1)], [y_b(1) y_b(1) y_b(1) y_b(1) y_b(1)], [z_b(1) z_b(1) z_b(2) z_b(2) z_b(1)], 'k-', 'LineWidth', 0.5);

scatter3(corners_x, corners_y, corners_z, 1, 'w', 'MarkerEdgeAlpha', 0);
for r = 1:2
    xyz = all_xyz{r};
    preds = all_preds_norm{r};
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, preds, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
end
colormap(custom_cmap); 
clim([-1, 1]); 
c3 = colorbar;
axis equal; grid off;
xlim(x_b); ylim(y_b); zlim(z_b);
set(gca, 'PlotBoxAspectRatio', [1 1 1], 'DataAspectRatio', [1 1 1]);
set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
set(gca, 'CameraPosition', cam_pos, 'CameraTarget', cam_tgt, 'CameraViewAngle', cam_va);
hold off;
saveas(fig3, fullfile(out_dir, 'Gradient_Hippocampus.png'));