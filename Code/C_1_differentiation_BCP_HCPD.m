clear; clc; 

Top = '/data/projects/punim2400';
addpath(genpath(fullfile(Top, 'tools')));

%% HCPD
amyg_idx = [2,10]; hipp_idx = [1,9];
cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
parc = cifti_label.cdata;
subcort_all_idx = find(ismember(parc, [hipp_idx, amyg_idx]));
n_vox_common = length(subcort_all_idx);
subcort_vals = parc(subcort_all_idx);
amyg_cols = ismember(subcort_vals, amyg_idx);
hipp_cols = ismember(subcort_vals, hipp_idx);

hcpd_input_dir = '/data/projects/punim2400/derivatives/HCPD/corr_maps';
hcpd_subs_master_table = readtable(fullfile(Top, 'HCPDfMRI', 'ndar_subject01.txt'));
hcpd_fd_table = readtable(fullfile(Top, 'derivatives', 'qc_metrics', 'mean_fd_allsubs.csv'));
hcpd_icv_table = readtable(fullfile(Top, 'derivatives', 'qc_metrics', 'icv_values.csv'));

hcpd_master_table = innerjoin(hcpd_subs_master_table, hcpd_fd_table, 'Keys', 'src_subject_id');
hcpd_master_table = innerjoin(hcpd_master_table, hcpd_icv_table, 'Keys', 'src_subject_id');

hcpd_files = dir(fullfile(hcpd_input_dir, 'CorrVals_AmygHipOnly_*.mat'));
n_subj_hcpd = length(hcpd_files);
delta_r_hcpd = nan(n_subj_hcpd, n_vox_common);
subj_ids_hcpd = cell(n_subj_hcpd, 1);

for s = 1:n_subj_hcpd
    D = load(fullfile(hcpd_input_dir, hcpd_files(s).name));
    subj_ids_hcpd{s} = D.subj_id; 
    
    c_amyg = zeros(1, n_vox_common); 
    c_hipp = zeros(1, n_vox_common);
    
    is_left = ismember(D.subcortical_labels, [9,10]); 
    is_right = ismember(D.subcortical_labels, [1,2]);
    
    c_amyg(is_left) = D.corr_left(is_left,2); c_hipp(is_left) = D.corr_left(is_left,1);
    c_amyg(is_right) = D.corr_right(is_right,2); c_hipp(is_right) = D.corr_right(is_right,1);
    
    delta_r_hcpd(s,:) = atanh(c_amyg) - atanh(c_hipp);
end

[~, master_idx] = ismember(subj_ids_hcpd, hcpd_master_table.src_subject_id);
aligned_cov_table = hcpd_master_table(master_idx(master_idx > 0), :);
aligned_delta_r = delta_r_hcpd(master_idx > 0, :);

fd_mask = aligned_cov_table.MeanFD <= 0.2;
hcpd_final_covariates = aligned_cov_table(fd_mask, :);
hcpd_final_delta_r = aligned_delta_r(fd_mask, :);

sex_numeric = nan(height(hcpd_final_covariates), 1);
for i = 1:height(hcpd_final_covariates)
    if strcmpi(hcpd_final_covariates.sex{i}, 'M')
        sex_numeric(i) = 1;
    elseif strcmpi(hcpd_final_covariates.sex{i}, 'F')
        sex_numeric(i) = 2;
    end
end
hcpd_final_covariates.sex = sex_numeric;

%% BCP
bcp_analysis_dir = fullfile(Top, 'derivatives', 'BCP');
bcp_input_dir = fullfile(bcp_analysis_dir, 'corr_maps');
bcp_master_table = readtable(fullfile(Top, 'BABY/image03/BCP', 'ndar_subject01.txt'));
bcp_master_table.src_subject_id = compose("sub-%06d", double(bcp_master_table.src_subject_id));

icv_table = readtable(fullfile(Top, 'derivatives', 'BCP_FC_Analysis', 'icv_values_BCP.csv'));
mean_fd_table = readtable(fullfile(Top, 'derivatives', 'BCP_FC_Analysis', 'mean_fd_run_BCP.csv')); 

icv_table.src_subject_id = string(icv_table.src_subject_id);
icv_table.session = string(icv_table.session);

mean_fd_table.src_subject_id = string(mean_fd_table.src_subject_id);
mean_fd_table.session = string(mean_fd_table.session);
mean_fd_table.run = string(mean_fd_table.run);

bcp_files = dir(fullfile(bcp_input_dir, 'CorrVals_AmygHipOnly_*.mat'));
n_scans_bcp = length(bcp_files);
delta_r_bcp = nan(n_scans_bcp, n_vox_common);

scan_info_subj = strings(n_scans_bcp, 1);
scan_info_sess = strings(n_scans_bcp, 1);
scan_info_run = strings(n_scans_bcp, 1); 
scan_ages_mo = nan(n_scans_bcp, 1);

for s = 1:n_scans_bcp
    D = load(fullfile(bcp_input_dir, bcp_files(s).name));
    raw_id = D.subj_id;
    if ischar(raw_id) || isstring(raw_id)
        id_num = str2double(raw_id);
    else
        id_num = double(raw_id);
    end
    
    scan_info_subj(s) = "sub-" + sprintf('%06d', id_num); 
    scan_info_sess(s) = string(D.ses_name);
    scan_info_run(s) = string(D.run);
    
    age_match = regexp(D.ses_name, '\d+', 'match');
    if ~isempty(age_match), scan_ages_mo(s) = str2double(age_match{1}); end
    
    c_amyg = zeros(1,n_vox_common); c_hipp = zeros(1,n_vox_common);
    is_left = ismember(D.subcortical_labels,[9,10]); is_right = ismember(D.subcortical_labels,[1,2]);
    c_amyg(is_left) = D.corr_left(is_left,2); c_hipp(is_left) = D.corr_left(is_left,1);
    c_amyg(is_right) = D.corr_right(is_right,2); c_hipp(is_right) = D.corr_right(is_right,1);
    
    delta_r_bcp(s,:) = atanh(c_amyg) - atanh(c_hipp);
end

scan_data_table = table(scan_info_subj, scan_info_sess, scan_info_run, scan_ages_mo, delta_r_bcp, ...
    'VariableNames', {'src_subject_id', 'session', 'run', 'interview_age', 'delta_r_raw'});

aligned_table = outerjoin(scan_data_table, bcp_master_table(:, {'src_subject_id','sex','race','site'}), ...
    'Keys','src_subject_id','MergeKeys',true,'Type','left');
aligned_table = outerjoin(aligned_table, icv_table, ...
    'Keys', {'src_subject_id', 'session'}, 'MergeKeys', true, 'Type', 'left');
aligned_table = outerjoin(aligned_table, mean_fd_table, ...
    'Keys', {'src_subject_id','session', 'run'}, 'MergeKeys', true, 'Type', 'left');

fd_threshold = 0.5;
drop_fd = aligned_table.MeanFD > fd_threshold | ismissing(aligned_table.MeanFD);
drop_nan_fc = any(isnan(aligned_table.delta_r_raw), 2);

essential_vars = {'interview_age', 'sex', 'site', 'MeanFD', 'ICV'};

temp_table = aligned_table(~drop_fd & ~drop_nan_fc, :);
keep_mask = true(height(temp_table), 1);

for i = 1:height(temp_table)
    for k = 1:length(essential_vars)
        val = temp_table.(essential_vars{k})(i);
        if ismissing(val)
            keep_mask(i) = false;
            break; 
        end
    end
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

[~, keep_indices_for_delta_r] = ismember(bcp_final_covariates(:,{'src_subject_id','session', 'run'}), ...
    scan_data_table(:,{'src_subject_id','session', 'run'}), 'rows');
bcp_final_delta_r = delta_r_bcp(keep_indices_for_delta_r(keep_indices_for_delta_r > 0),:);

bcp_final_covariates.delta_r_raw = [];

sex_numeric = nan(height(bcp_final_covariates), 1);
for i = 1:height(bcp_final_covariates)
    val = bcp_final_covariates.sex(i); 
    if strcmpi(val, "M")
        sex_numeric(i) = 1;
    elseif strcmpi(val, "F")
        sex_numeric(i) = 2;
    end
end
bcp_final_covariates.sex = sex_numeric;

%% Harmonize
cols_to_keep = {'src_subject_id', 'interview_age', 'sex', 'site', 'MeanFD', 'ICV'};
master_covariates = [hcpd_final_covariates(:, cols_to_keep); ...
                     bcp_final_covariates(:, cols_to_keep)];
master_covariates.interview_age = master_covariates.interview_age / 12;

[master_covariates, removed_rows_idx] = rmmissing(master_covariates);
master_delta_r = [hcpd_final_delta_r; bcp_final_delta_r];
master_delta_r(removed_rows_idx,:) = [];
master_delta_r = master_delta_r(:,1:1743);

dat_to_harmonize = master_delta_r';
[~, ~, batch] = unique(master_covariates.site);

age_z = zscore(master_covariates.interview_age);
sex_dummy = dummyvar(master_covariates.sex);
meanfd_z = zscore(master_covariates.MeanFD);
icv_z = zscore(master_covariates.ICV);

mod = [age_z, sex_dummy(:, 1:end-1), meanfd_z, icv_z];

harmonized_dat = combat(dat_to_harmonize, batch, mod, 1);
harmonized_delta_r_voxels = harmonized_dat';

output_raw = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_voxel_BCPHCPD.csv');
writematrix(harmonized_dat', output_raw);

deltar_a_harmonized = mean(harmonized_delta_r_voxels(:, amyg_cols), 2, 'omitnan');
deltar_h_harmonized = mean(harmonized_delta_r_voxels(:, hipp_cols), 2, 'omitnan');

harmonized_master_table = master_covariates;
harmonized_master_table.deltar_a = deltar_a_harmonized;
harmonized_master_table.deltar_h = deltar_h_harmonized;

writetable(harmonized_master_table, fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_table_BCPHCPD.csv'));