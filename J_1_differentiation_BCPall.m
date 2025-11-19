clear; clc; 
Top = '/data/project';
addpath(genpath(fullfile(Top, 'tools')));


amyg_idx = [2,10]; hipp_idx = [1,9];
cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
parc = cifti_label.cdata;
amyg_hip_vox_idx = find(ismember(parc, [hipp_idx, amyg_idx]));
n_vox_common = length(amyg_hip_vox_idx);

bcp_analysis_dir = fullfile(Top, 'derivatives', 'BCP');
bcp_input_dir = fullfile(bcp_analysis_dir, 'corr_maps');

bcp_master_table = readtable(fullfile(Top, 'BABY/image03/BCP', 'ndar_subject01.txt'));
if isnumeric(bcp_master_table.src_subject_id), bcp_master_table.src_subject_id = arrayfun(@(x) sprintf('sub-%06d', x), bcp_master_table.src_subject_id, 'UniformOutput', false); end

icv_table = readtable(fullfile(bcp_analysis_dir, 'icv_values_BCP.csv'));
mean_fd_table = readtable(fullfile(bcp_analysis_dir, 'mean_fd_run_BCP.csv')); 
if isnumeric(mean_fd_table.src_subject_id), mean_fd_table.src_subject_id = arrayfun(@(x) sprintf('sub-%06d', x), mean_fd_table.src_subject_id, 'UniformOutput', false); end
if ~isstring(mean_fd_table.run)
    mean_fd_table.run = string(mean_fd_table.run);
end

bcp_files = dir(fullfile(bcp_input_dir, 'CorrVals_AmygHipOnly_*.mat'));
n_scans_bcp = length(bcp_files);
delta_r_bcp = nan(n_scans_bcp, n_vox_common);
scan_info_subj = cell(n_scans_bcp, 1);
scan_info_sess = cell(n_scans_bcp, 1);
scan_info_run = strings(n_scans_bcp, 1); 
scan_ages_mo = nan(n_scans_bcp, 1);

for s = 1:n_scans_bcp
    D = load(fullfile(bcp_input_dir, bcp_files(s).name));
    
    scan_info_subj{s} = ['sub-' D.subj_id]; 
    scan_info_sess{s} = D.ses_name;
    scan_info_run(s) = string(D.run);
    
    age_match = regexp(D.ses_name, '\d+', 'match');
    if ~isempty(age_match), scan_ages_mo(s) = str2double(age_match{1}); end
    
    c_amyg=zeros(1,n_vox_common); c_hipp=zeros(1,n_vox_common);
    is_left=ismember(D.subcortical_labels,[9,10]); is_right=ismember(D.subcortical_labels,[1,2]);
    c_amyg(is_left)=D.corr_left(is_left,2); c_hipp(is_left)=D.corr_left(is_left,1);
    c_amyg(is_right)=D.corr_right(is_right,2); c_hipp(is_right)=D.corr_right(is_right,1);
    
    delta_r_bcp(s,:) = atanh(c_amyg) - atanh(c_hipp);
end

scan_data_table = table(scan_info_subj, scan_info_sess, scan_info_run, scan_ages_mo, delta_r_bcp, ...
    'VariableNames', {'src_subject_id', 'session', 'run', 'interview_age', 'delta_r_raw'});

aligned_table = outerjoin(scan_data_table, bcp_master_table(:, {'src_subject_id','sex','race','site'}), 'Keys','src_subject_id','MergeKeys',true,'Type','left');
aligned_table = outerjoin(aligned_table, icv_table, 'Keys', {'src_subject_id', 'session'}, 'MergeKeys', true, 'Type', 'left');
aligned_table = outerjoin(aligned_table, mean_fd_table, 'Keys', {'src_subject_id','session', 'run'}, 'MergeKeys', true, 'Type', 'left');

fd_threshold = 0.5;
drop_fd = aligned_table.MeanFD > fd_threshold | ismissing(aligned_table.MeanFD);
drop_nan_fc = any(isnan(aligned_table.delta_r_raw), 2);
temp_table = aligned_table(~drop_fd & ~drop_nan_fc, :);

essential_vars = {'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV'};
keep_mask = true(height(temp_table), 1);
for k = 1:numel(essential_vars)
    v = temp_table.(essential_vars{k});
    if iscell(v) || iscategorical(v)
        if iscell(v), keep_mask = keep_mask & ~cellfun(@(x) isempty(x) || (ischar(x) && isempty(strtrim(x))), v);
        else, keep_mask = keep_mask & ~isundefined(v); end
    else, keep_mask = keep_mask & ~ismissing(v); end
end
bcp_covariates_filtered = temp_table(keep_mask, :);

sorted_table = sortrows(bcp_covariates_filtered, {'src_subject_id', 'MeanFD'});
[~, best_run_indices] = unique(sorted_table(:, {'src_subject_id', 'session'}), 'rows', 'first');
best_sessions_table = sorted_table(best_run_indices, :);

sorted_by_age_table = sortrows(best_sessions_table, ...
                               {'src_subject_id', 'interview_age'}, ...
                               {'ascend', 'descend'});
[~, oldest_session_indices] = unique(sorted_by_age_table.src_subject_id, 'first');
bcp_final_covariates = sorted_by_age_table(oldest_session_indices, :);

[~, keep_indices_for_delta_r] = ismember(bcp_final_covariates(:,{'src_subject_id','session', 'run'}), scan_data_table(:,{'src_subject_id','session', 'run'}), 'rows');
bcp_final_delta_r = delta_r_bcp(keep_indices_for_delta_r(keep_indices_for_delta_r > 0),:);
bcp_final_covariates.delta_r_raw = [];

sex_numeric = nan(height(bcp_final_covariates), 1);
sex_numeric(strcmpi(bcp_final_covariates.sex, 'M')) = 1;
sex_numeric(strcmpi(bcp_final_covariates.sex, 'F')) = 2;
bcp_final_covariates.sex = sex_numeric;
[~,~,bcp_final_covariates.race] = unique(bcp_final_covariates.race);

cols_to_keep = {'src_subject_id', 'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV'};
master_covariates = bcp_final_covariates(:, cols_to_keep);
master_covariates.interview_age = master_covariates.interview_age / 12;

[master_covariates, removed_rows_idx] = rmmissing(master_covariates);
master_delta_r = bcp_final_delta_r;
master_delta_r(removed_rows_idx,:) = [];

dat_to_harmonize = master_delta_r';
[~, ~, batch] = unique(master_covariates.site);

age_z = zscore(master_covariates.interview_age);
sex_dummy = dummyvar(master_covariates.sex);
race_dummy = dummyvar(master_covariates.race);
meanfd_z = zscore(master_covariates.MeanFD);
icv_z = zscore(master_covariates.ICV);

mod = [age_z, sex_dummy(:, 1:end-1), race_dummy(:, 1:end-1), meanfd_z, icv_z];

harmonized_dat = combat(dat_to_harmonize, batch, mod, 1);
harmonized_delta_r_voxels = harmonized_dat';

output_raw = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_voxel_BCPall.csv');
writematrix(harmonized_dat', output_raw);

deltar_a_harmonized = mean(harmonized_delta_r_voxels(:, amyg_cols), 2, 'omitnan');
deltar_h_harmonized = mean(harmonized_delta_r_voxels(:, hipp_cols), 2, 'omitnan');

harmonized_master_table = master_covariates;
harmonized_master_table.deltar_a = deltar_a_harmonized;
harmonized_master_table.deltar_h = deltar_h_harmonized;

output_filename = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_table_BCPall.csv');
writetable(harmonized_master_table, output_filename);
