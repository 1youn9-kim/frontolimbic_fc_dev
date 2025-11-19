clear; clc;


Top = '/data/project';
Path = fullfile(Top, 'PNC/Data'); 
addpath(genpath(fullfile(Top, 'tools')));
output_dir = fullfile(Top, 'derivatives/binned_fc_maps_amyg_PNC');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

subcortical_r = [1, 2]; subcortical_l = [9, 10];
medial = [42,57,59,73,74,75,76,77,78,79,80,81,85,86,88,104,106,109,180,181,182,195,196,222,237,239,253,254,255,256,257,258,259,260,261,265,266,268,284,286,289,360,361,362,375,376];
lateral = [26,27,28,60,84,87,93,94,96,97,98,100,101,102,103,105,107,112,113,186,187,206,207,208,240,264,267,273,274,277,278,280,281,282,283,285,287,292,293,366,367,82,83,89,90,91,92,95,99,108,110,114,262,263,269,270,271,272,275,276,279,288,290,294];
ROI_frontal_cortex = [medial, lateral];
atlas_path = fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
cifti_label = cifti_read(atlas_path); parc = cifti_label.cdata;
ctx_vertex_idx = find(ismember(parc, ROI_frontal_cortex));
medial_idx = find(ismember(ctx_vertex_idx, find(ismember(parc,medial))));
lateral_idx = find(ismember(ctx_vertex_idx, find(ismember(parc,lateral))));
subcort_all_idx = find(parc >= 1 & parc <= 16);
amyg_mask_sub = ismember(parc(subcort_all_idx), [2, 10]);

HCPref = load(fullfile(Top, 'derivatives/FC_reference_HCP/FC_reference_HCP_scan1.mat'));
ref = HCPref.FC_mean; 
ref = ref(:, ctx_vertex_idx);
ref_parc_subcort = parc(subcort_all_idx);
ref_amyg_l = mean(ref(ismember(ref_parc_subcort, 10), :), 1);
ref_amyg_r = mean(ref(ismember(ref_parc_subcort, 2), :), 1);
ref_vec = (ref_amyg_l + ref_amyg_r) / 2;

subs_master_table = readtable(fullfile(Top, 'PNC_subjectlist.csv'));
subs_master_table = renamevars(subs_master_table, 'bids', 'src_subject_id');
if isnumeric(subs_master_table.src_subject_id), subs_master_table.src_subject_id = cellstr(num2str(subs_master_table.src_subject_id)); end
gamlss_table_pnc = readtable(fullfile(Top, 'derivatives/gamlss_scores/deltar_gamlss_zscores_PNC.csv'));
if isnumeric(gamlss_table_pnc.src_subject_id), gamlss_table_pnc.src_subject_id = cellstr(num2str(gamlss_table_pnc.src_subject_id)); end
zscore_col_name = 'deltar_a_zscore';
[~, ia, ib] = intersect(subs_master_table.src_subject_id, gamlss_table_pnc.src_subject_id);
subs = [subs_master_table(ia,:), gamlss_table_pnc(ib, {zscore_col_name})];

raw_fc_data_r = nan(height(subs), length(ctx_vertex_idx));
valid_subj_mask = true(height(subs), 1);

for s = 1:height(subs)
    subj_id = subs.src_subject_id{s};
    ts_path = fullfile(Path, subj_id, 'func', sprintf('%s_task-rest_space-fsLR_den-91k_bold.dtseries.nii', subj_id));
    conf_path = fullfile(Path, subj_id, 'func', sprintf('%s_task-rest_desc-confounds_timeseries.tsv', subj_id));
    if ~isfile(ts_path) || ~isfile(conf_path)
        valid_subj_mask(s) = false; continue;
    end
    
    bold_data_full = cifti_read(ts_path).cdata';
    conf_table = readtable(conf_path, 'FileType', 'text', 'Delimiter', '\t');
    fd = conf_table.framewise_displacement;
    scrub_idx = find(fd > 0.5);
    frames_to_censor = unique([scrub_idx-1; scrub_idx; scrub_idx+1; scrub_idx+2]);
    frames_to_censor(frames_to_censor < 1 | frames_to_censor > size(bold_data_full, 1)) = [];
    
    bold_data = bold_data_full;
    if ~isempty(frames_to_censor)
        good_frames = setdiff(1:size(bold_data_full, 1), frames_to_censor);
        for v = 1:size(bold_data_full, 2)
            bold_data(frames_to_censor, v) = interp1(good_frames, bold_data_full(good_frames, v), frames_to_censor, 'linear', 'extrap');
        end
    end
    bold_data = bold_data(5:end, :);
    conf_table = conf_table(5:end, :);
    
    motion_base = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
    motion_deriv = cellfun(@(x) [x, '_derivative1'], motion_base, 'UniformOutput', false);
    motion_24p = [motion_base, motion_deriv, cellfun(@(x) [x, '_power2'], motion_base, 'UniformOutput', false), cellfun(@(x) [x, '_power2'], motion_deriv, 'UniformOutput', false)];
    dct_cols = conf_table.Properties.VariableNames(startsWith(conf_table.Properties.VariableNames, 'cosine'));
    nuisance = table2array(conf_table(:, [motion_24p, {'white_matter', 'csf'}, dct_cols]));
    nuisance(isnan(nuisance)) = 0;
    
    residuals = bold_data - [ones(size(nuisance, 1), 1), nuisance] * ([ones(size(nuisance, 1), 1), nuisance] \ bold_data);
    A = normalize(residuals', 2);
    FC_vox_profile = corr(A(subcort_all_idx, :)', A(ctx_vertex_idx, :)');
    raw_fc_data_r(s, :) = mean(FC_vox_profile(amyg_mask_sub, :), 1);
end
subs = subs(valid_subj_mask, :);
raw_fc_data_r = raw_fc_data_r(valid_subj_mask, :);
raw_fc_data_r(raw_fc_data_r > 1) = 1; raw_fc_data_r(raw_fc_data_r < -1) = -1;
raw_fc_data_z = atanh(raw_fc_data_r);


mean_row = mean(raw_fc_data_z, 2, 'omitnan');
std_row  = std(raw_fc_data_z, 0, 2, 'omitnan'); std_row(std_row == 0) = 1; 
harmonized_data_z = (raw_fc_data_z - mean_row) ./ std_row;

%% Age Bin
age_bins = 8:21;
template = cifti_label; template.cdata = nan(size(cifti_label.cdata,1), 1); template.diminfo{2} = cifti_diminfo_make_scalars(1);
corrmat_m = zeros(length(age_bins),3); corrmat_l = zeros(length(age_bins),3); corrmat_c = zeros(length(age_bins),3);
corrmat_m_all = zeros(length(age_bins),1); corrmat_l_all = zeros(length(age_bins),1);

for i = 1:length(age_bins)
    idx = find((subs.age >= age_bins(i)) & (subs.age < age_bins(i) + 1));
    if isempty(idx), continue; end
    
    final_map_all = tanh(mean(harmonized_data_z(idx, :), 1, 'omitnan'))';
    out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map_all;
    cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_AgeBin%02d_all.dscalar.nii', age_bins(i))));
    corrmat_m_all(i) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map_all(medial_idx), 'rows','complete');
    corrmat_l_all(i) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map_all(lateral_idx), 'rows','complete');
    
    z_scr = subs.(zscore_col_name)(idx);
    if length(z_scr) < 3, continue; end
    beh_tertiles = quantile(z_scr, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scr < beh_tertiles(1)), ...
                   'mid', idx(z_scr >= beh_tertiles(1) & z_scr < beh_tertiles(2)), ...
                   'high', idx(z_scr >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labs = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    
    for g = 1:3
        s_idx = g_idx.(g_keys{g});
        if isempty(s_idx), continue; end
        final_map = tanh(mean(harmonized_data_z(s_idx, :), 1, 'omitnan'))';
        out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map;
        cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_AgeBin%02d_%s.dscalar.nii', age_bins(i), g_labs{g})));
        corrmat_m(i,g) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map(medial_idx), 'rows','complete');
        corrmat_l(i,g) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map(lateral_idx), 'rows','complete');
        corrmat_c(i,g) = corr(ref_AvH_Tmap_frontal, final_map, 'rows','complete');
    end
end
writematrix(corrmat_m, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBinYearly_CorrWRef_medial.csv'));
writematrix(corrmat_l, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBinYearly_CorrWRef_lateral.csv'));
writematrix(corrmat_c, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBinYearly_CorrWRef.csv'));
writematrix(corrmat_m_all, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBinYearly_CorrWRef_medial_all.csv'));
writematrix(corrmat_l_all, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBinYearly_CorrWRef_lateral_all.csv'));

%% Tertiles
fprintf('Binning (Age Tertiles)...\n');
at = quantile(subs.age, [1/3, 2/3]);
masks = {subs.age < at(1), subs.age >= at(1) & subs.age < at(2), subs.age >= at(2)};
labs = {'Age_Tertile1', 'Age_Tertile2', 'Age_Tertile3'};
corrmat_m_3 = zeros(3,3); corrmat_l_3 = zeros(3,3); corrmat_c_3 = zeros(3,3);
corrmat_m_all_3 = zeros(3,1); corrmat_l_all_3 = zeros(3,1);

for i = 1:3
    idx = find(masks{i}); if isempty(idx), continue; end
    final_map_all = tanh(mean(harmonized_data_z(idx, :), 1, 'omitnan'))';
    out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map_all;
    cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_Denoised_%s_all.dscalar.nii', labs{i})));
    corrmat_m_all_3(i) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map_all(medial_idx), 'rows','complete');
    corrmat_l_all_3(i) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map_all(lateral_idx), 'rows','complete');
    
    z_scr = subs.(zscore_col_name)(idx);
    beh_tertiles = quantile(z_scr, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scr < beh_tertiles(1)), ...
                   'mid', idx(z_scr >= beh_tertiles(1) & z_scr < beh_tertiles(2)), ...
                   'high', idx(z_scr >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labs = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    
    for g = 1:3
        s_idx = g_idx.(g_keys{g}); if isempty(s_idx), continue; end
        final_map = tanh(mean(harmonized_data_z(s_idx, :), 1, 'omitnan'))';
        out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map;
        cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_Denoised_%s_%s.dscalar.nii', labs{i}, g_labs{g})));
        corrmat_m_3(i,g) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map(medial_idx), 'rows','complete');
        corrmat_l_3(i,g) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map(lateral_idx), 'rows','complete');
        corrmat_c_3(i,g) = corr(ref_AvH_Tmap_frontal, final_map, 'rows','complete');
    end
end
writematrix(corrmat_m_3, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin3Split_CorrWRef_medial.csv'));
writematrix(corrmat_l_3, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin3Split_CorrWRef_lateral.csv'));
writematrix(corrmat_c_3, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin3Split_CorrWRef.csv'));
writematrix(corrmat_m_all_3, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin3Split_CorrWRef_medial_all.csv'));
writematrix(corrmat_l_all_3, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin3Split_CorrWRef_lateral_all.csv'));

%% Quartiles
fprintf('Binning (Age Quartiles)...\n');
aq = quantile(subs.age, [0.25, 0.5, 0.75]);
masks = {subs.age < aq(1), subs.age >= aq(1) & subs.age < aq(2), subs.age >= aq(2) & subs.age < aq(3), subs.age >= aq(3)};
labs = {'Age_Quartile1', 'Age_Quartile2', 'Age_Quartile3', 'Age_Quartile4'};
corrmat_m_4 = zeros(4,3); corrmat_l_4 = zeros(4,3); corrmat_c_4 = zeros(4,3);
corrmat_m_all_4 = zeros(4,1); corrmat_l_all_4 = zeros(4,1);

for i = 1:4
    idx = find(masks{i}); if isempty(idx), continue; end
    final_map_all = tanh(mean(harmonized_data_z(idx, :), 1, 'omitnan'))';
    out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map_all;
    cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_Denoised_%s_all.dscalar.nii', labs{i})));
    corrmat_m_all_4(i) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map_all(medial_idx), 'rows','complete');
    corrmat_l_all_4(i) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map_all(lateral_idx), 'rows','complete');
    
    z_scr = subs.(zscore_col_name)(idx);
    beh_tertiles = quantile(z_scr, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scr < beh_tertiles(1)), ...
                   'mid', idx(z_scr >= beh_tertiles(1) & z_scr < beh_tertiles(2)), ...
                   'high', idx(z_scr >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labs = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    
    for g = 1:3
        s_idx = g_idx.(g_keys{g}); if isempty(s_idx), continue; end
        final_map = tanh(mean(harmonized_data_z(s_idx, :), 1, 'omitnan'))';
        out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map;
        cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_Denoised_%s_%s.dscalar.nii', labs{i}, g_labs{g})));
        corrmat_m_4(i,g) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map(medial_idx), 'rows','complete');
        corrmat_l_4(i,g) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map(lateral_idx), 'rows','complete');
        corrmat_c_4(i,g) = corr(ref_AvH_Tmap_frontal, final_map, 'rows','complete');
    end
end
writematrix(corrmat_m_4, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin4Split_CorrWRef_medial.csv'));
writematrix(corrmat_l_4, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin4Split_CorrWRef_lateral.csv'));
writematrix(corrmat_c_4, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin4Split_CorrWRef.csv'));
writematrix(corrmat_m_all_4, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin4Split_CorrWRef_medial_all.csv'));
writematrix(corrmat_l_all_4, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin4Split_CorrWRef_lateral_all.csv'));

%% Median Split
fprintf('Binning (Median Split)...\n');
am = quantile(subs.age, 0.5);
masks = {subs.age < am, subs.age >= am};
labs = {'Age_MedianSplit_Low', 'Age_MedianSplit_High'};
corrmat_m_2 = zeros(2,3); corrmat_l_2 = zeros(2,3); corrmat_c_2 = zeros(2,3);
corrmat_m_all_2 = zeros(2,1); corrmat_l_all_2 = zeros(2,1);

for i = 1:2
    idx = find(masks{i}); if isempty(idx), continue; end
    final_map_all = tanh(mean(harmonized_data_z(idx, :), 1, 'omitnan'))';
    out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map_all;
    cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_Denoised_%s_all.dscalar.nii', labs{i})));
    corrmat_m_all_2(i) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map_all(medial_idx), 'rows','complete');
    corrmat_l_all_2(i) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map_all(lateral_idx), 'rows','complete');
    
    z_scr = subs.(zscore_col_name)(idx);
    beh_tertiles = quantile(z_scr, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scr < beh_tertiles(1)), ...
                   'mid', idx(z_scr >= beh_tertiles(1) & z_scr < beh_tertiles(2)), ...
                   'high', idx(z_scr >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labs = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    
    for g = 1:3
        s_idx = g_idx.(g_keys{g}); if isempty(s_idx), continue; end
        final_map = tanh(mean(harmonized_data_z(s_idx, :), 1, 'omitnan'))';
        out_cifti = template; out_cifti.cdata(ctx_vertex_idx) = final_map;
        cifti_write(out_cifti, fullfile(output_dir, sprintf('AmygAvgFC_PNC_Denoised_%s_%s.dscalar.nii', labs{i}, g_labs{g})));
        corrmat_m_2(i,g) = corr(ref_AvH_Tmap_frontal(medial_idx), final_map(medial_idx), 'rows','complete');
        corrmat_l_2(i,g) = corr(ref_AvH_Tmap_frontal(lateral_idx), final_map(lateral_idx), 'rows','complete');
        corrmat_c_2(i,g) = corr(ref_AvH_Tmap_frontal, final_map, 'rows','complete');
    end
end
writematrix(corrmat_m_2, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin2Split_CorrWRef_medial.csv'));
writematrix(corrmat_l_2, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin2Split_CorrWRef_lateral.csv'));
writematrix(corrmat_c_2, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin2Split_CorrWRef.csv'));
writematrix(corrmat_m_all_2, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin2Split_CorrWRef_medial_all.csv'));
writematrix(corrmat_l_all_2, fullfile(output_dir, 'AmygAvgFC_PNC_Denoised_AgeBin2Split_CorrWRef_lateral_all.csv'));
fprintf('Done.\n');