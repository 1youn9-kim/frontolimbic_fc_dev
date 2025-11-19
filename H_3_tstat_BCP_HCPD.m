clear; clc; 

Top = '/data/projects/punim2400';
addpath(genpath(fullfile(Top, 'tools')));

amyg_idx = [2,10]; hipp_idx = [1,9];
cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
parc = cifti_label.cdata;
amyg_hip_vox_idx = find(ismember(parc, [hipp_idx, amyg_idx]));
n_vox_common = length(amyg_hip_vox_idx);

amyg_cols = ismember(parc(amyg_hip_vox_idx), amyg_idx); 
hipp_cols = ismember(parc(amyg_hip_vox_idx), hipp_idx);  
alpha = 0.05;

%% HCPD
hcpd_input_dir = fullfile(Top, 'derivatives/HCPD/tstat_maps_final');
hcpd_subs_master_table = readtable(fullfile(Top, 'HCPDfMRI/ndar_subject01.txt'));
hcpd_fd_table = readtable(fullfile(Top, 'derivatives/qc_metrics/mean_fd_allsubs.csv'));
hcpd_icv_table = readtable(fullfile(Top, 'derivatives/qc_metrics/icv_values.csv'));

hcpd_master_table = innerjoin(hcpd_subs_master_table, hcpd_fd_table, 'Keys', 'src_subject_id');
hcpd_master_table = innerjoin(hcpd_master_table, hcpd_icv_table, 'Keys', 'src_subject_id');

hcpd_files = dir(fullfile(hcpd_input_dir, 'Tstat_AmygHipCompare_*.mat'));
n_subj_hcpd = length(hcpd_files);
delta_r_hcpd = nan(n_subj_hcpd, n_vox_common);
direg_scores_hcpd = nan(n_subj_hcpd, 2);
subj_ids_hcpd = cell(n_subj_hcpd, 1);

for s = 1:n_subj_hcpd
    D = load(fullfile(hcpd_input_dir, hcpd_files(s).name));
    if ~isequal(D.voxel_indices, amyg_hip_vox_idx), error('Voxel space mismatch!'); end
    subj_ids_hcpd{s} = D.subj_id; 
    
    t = D.t_stat_map;
    p = D.p_vals;
    delta_r_hcpd(s,:) = t';
    
    mean_good_amyg = mean(t(p<alpha & t>0 & amyg_cols));
    mean_bad_amyg  = mean(abs(t(p<alpha & t<0 & amyg_cols)));
    if isnan(mean_good_amyg), mean_good_amyg = 0; end
    if isnan(mean_bad_amyg), mean_bad_amyg = 0; end
    
    mean_good_hipp = mean(abs(t(p<alpha & t<0 & hipp_cols)));
    mean_bad_hipp  = mean(t(p<alpha & t>0 & hipp_cols));
    if isnan(mean_good_hipp), mean_good_hipp = 0; end
    if isnan(mean_bad_hipp), mean_bad_hipp = 0; end
    
    direg_scores_hcpd(s, :) = [(mean_good_amyg - mean_bad_amyg), (mean_good_hipp - mean_bad_hipp)];
end

[~, master_idx] = ismember(subj_ids_hcpd, hcpd_master_table.src_subject_id);
aligned_cov_table = hcpd_master_table(master_idx(master_idx > 0), :);
aligned_delta_r = delta_r_hcpd(master_idx > 0, :);
aligned_direg = direg_scores_hcpd(master_idx > 0, :);

fd_mask = aligned_cov_table.MeanFD <= 0.2;
hcpd_final_covariates = aligned_cov_table(fd_mask, :);
hcpd_final_delta_r = aligned_delta_r(fd_mask, :);
hcpd_final_direg = aligned_direg(fd_mask, :);

sex_numeric = nan(height(hcpd_final_covariates), 1);
for i = 1:height(hcpd_final_covariates)
    if strcmpi(hcpd_final_covariates.sex{i}, 'M')
        sex_numeric(i) = 1;
    elseif strcmpi(hcpd_final_covariates.sex{i}, 'F')
        sex_numeric(i) = 2;
    end
end
hcpd_final_covariates.sex = sex_numeric;
[~,~,hcpd_final_covariates.race] = unique(hcpd_final_covariates.race);

%% BCP
bcp_analysis_dir = fullfile(Top, 'derivatives', 'BCP', 'tstat_maps_final');
bcp_master_table = readtable(fullfile(Top, 'BABY/image03/BCP', 'ndar_subject01.txt'));
if isnumeric(bcp_master_table.src_subject_id)
    bcp_master_table.src_subject_id = arrayfun(@(x) sprintf('sub-%06d', x), bcp_master_table.src_subject_id, 'UniformOutput', false);
end

icv_table = readtable(fullfile(Top, 'derivatives', 'BCP', 'icv_values_BCP.csv'));
mean_fd_table = readtable(fullfile(Top, 'derivatives', 'BCP', 'mean_fd_run_BCP.csv')); 
if isnumeric(mean_fd_table.src_subject_id)
    mean_fd_table.src_subject_id = arrayfun(@(x) sprintf('sub-%06d', x), mean_fd_table.src_subject_id, 'UniformOutput', false);
end
mean_fd_table.run = string(mean_fd_table.run);

bcp_files = dir(fullfile(bcp_analysis_dir, 'Tstat_AmygHipCompare_*.mat'));
n_scans_bcp = length(bcp_files);
delta_r_bcp = nan(n_scans_bcp, n_vox_common);
direg_scores_bcp = nan(n_scans_bcp, 2);
scan_info_subj = cell(n_scans_bcp, 1);
scan_info_sess = cell(n_scans_bcp, 1);
scan_info_run = strings(n_scans_bcp, 1); 
scan_ages_mo = nan(n_scans_bcp, 1);

for s = 1:n_scans_bcp
    D = load(fullfile(bcp_analysis_dir, bcp_files(s).name));
    if ~isequal(D.voxel_indices, amyg_hip_vox_idx), error('Voxel space mismatch!'); end
    
    scan_info_subj{s} = ['sub-' D.subj_id]; 
    scan_info_sess{s} = D.ses_name;
    scan_info_run(s) = string(D.run);
    
    age_match = regexp(D.ses_name, '\d+', 'match');
    if ~isempty(age_match), scan_ages_mo(s) = str2double(age_match{1}); end
    
    t = D.t_stat_map;
    p = D.p_vals;
    delta_r_bcp(s,:) = t';
    
    mean_good_amyg = mean(t(p<alpha & t>0 & amyg_cols));
    mean_bad_amyg  = mean(abs(t(p<alpha & t<0 & amyg_cols)));
    if isnan(mean_good_amyg), mean_good_amyg = 0; end
    if isnan(mean_bad_amyg), mean_bad_amyg = 0; end
    
    mean_good_hipp = mean(abs(t(p<alpha & t<0 & hipp_cols)));
    mean_bad_hipp  = mean(t(p<alpha & t>0 & hipp_cols));
    if isnan(mean_good_hipp), mean_good_hipp = 0; end
    if isnan(mean_bad_hipp), mean_bad_hipp = 0; end
    
    direg_scores_bcp(s, :) = [(mean_good_amyg - mean_bad_amyg), (mean_good_hipp - mean_bad_hipp)];
end

scan_data_table = table(scan_info_subj, scan_info_sess, scan_info_run, scan_ages_mo, delta_r_bcp, direg_scores_bcp, ...
    'VariableNames', {'src_subject_id', 'session', 'run', 'interview_age', 'delta_r_raw', 'direg_raw'});

aligned_table = outerjoin(scan_data_table, bcp_master_table(:, {'src_subject_id','sex','race','site'}), ...
    'Keys','src_subject_id','MergeKeys',true,'Type','left');
aligned_table = outerjoin(aligned_table, icv_table, ...
    'Keys', {'src_subject_id', 'session'}, 'MergeKeys', true, 'Type', 'left');
aligned_table = outerjoin(aligned_table, mean_fd_table, ...
    'Keys', {'src_subject_id','session', 'run'}, 'MergeKeys', true, 'Type', 'left');

fd_threshold = 0.5;
drop_fd = aligned_table.MeanFD > fd_threshold | ismissing(aligned_table.MeanFD);
drop_nan_fc = any(isnan(aligned_table.delta_r_raw), 2);

essential_vars = {'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV'};
temp_table = aligned_table(~drop_fd & ~drop_nan_fc, :);
keep_mask = true(height(temp_table), 1);

for i = 1:height(temp_table)
    for k = 1:length(essential_vars)
        val = temp_table.(essential_vars{k})(i);
        if ismissing(val), keep_mask(i) = false; break; end
    end
end

bcp_covariates_filtered = temp_table(keep_mask, :);

sorted_table = sortrows(bcp_covariates_filtered, {'src_subject_id', 'MeanFD'});
[~, best_run_indices] = unique(sorted_table(:, {'src_subject_id', 'session'}), 'rows', 'first');
best_sessions_table = sorted_table(best_run_indices, :);

sorted_by_age_table = sortrows(best_sessions_table, {'src_subject_id', 'interview_age'}, {'ascend', 'descend'});
[~, oldest_session_indices] = unique(sorted_by_age_table.src_subject_id, 'first');
bcp_final_covariates = sorted_by_age_table(oldest_session_indices, :);

bcp_final_delta_r = bcp_final_covariates.delta_r_raw;
bcp_final_direg = bcp_final_covariates.direg_raw;
bcp_final_covariates.delta_r_raw = [];
bcp_final_covariates.direg_raw = [];

sex_numeric = nan(height(bcp_final_covariates), 1);
for i = 1:height(bcp_final_covariates)
    val = bcp_final_covariates.sex{i};
    if strcmpi(val, 'M'), sex_numeric(i) = 1;
    elseif strcmpi(val, 'F'), sex_numeric(i) = 2; end
end
bcp_final_covariates.sex = sex_numeric;
[~,~,bcp_final_covariates.race] = unique(bcp_final_covariates.race);

%% Harmonize
cols_to_keep = {'src_subject_id', 'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV'};
master_covariates = [hcpd_final_covariates(:, cols_to_keep); bcp_final_covariates(:, cols_to_keep)];
master_covariates.interview_age = master_covariates.interview_age / 12;

master_delta_r = [hcpd_final_delta_r; bcp_final_delta_r];
master_direg = [hcpd_final_direg; bcp_final_direg];

[master_covariates, removed_rows_idx] = rmmissing(master_covariates);
master_delta_r(removed_rows_idx,:) = [];
master_direg(removed_rows_idx,:) = [];

[~, ~, batch] = unique(master_covariates.site);
age_z = zscore(master_covariates.interview_age);
sex_dummy = dummyvar(master_covariates.sex);
race_dummy = dummyvar(master_covariates.race);
meanfd_z = zscore(master_covariates.MeanFD);
icv_z = zscore(master_covariates.ICV);

mod = [age_z, sex_dummy(:, 1:end-1), race_dummy(:, 1:end-1), meanfd_z, icv_z];

harmonized_dat_voxels = combat(master_delta_r', batch, mod, 1)';
harmonized_dat_scores = combat(master_direg', batch, mod, 1)';

output_raw = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_DIreg_voxel_BCPHCPD.csv');
writematrix(harmonized_dat_voxels, output_raw);

harmonized_master_table = master_covariates;
harmonized_master_table.deltar_a = harmonized_dat_scores(:, 1);
harmonized_master_table.deltar_h = harmonized_dat_scores(:, 2);

writetable(harmonized_master_table, fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_DIreg_table_BCPHCPD.csv'));