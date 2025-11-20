clear; clc; 


Top = '/data/project';
addpath(genpath(fullfile(Top, 'tools')));


amyg_idx = [2,10]; hipp_idx = [1,9];
cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
parc = cifti_label.cdata;
amyg_hip_vox_idx = find(ismember(parc, [hipp_idx, amyg_idx]));
n_vox_common = length(amyg_hip_vox_idx);

hcpd_input_dir = '/data/projects/punim2400/derivatives/HCPD/corr_maps';

hcpd_subs_master_table = readtable(fullfile(Top, 'HCPDfMRI/ndar_subject01.txt'));
hcpd_fd_table = readtable(fullfile(Top, 'derivatives/qc_metrics/mean_fd_allsubs.csv'));
hcpd_icv_table = readtable(fullfile(Top, 'derivatives/qc_metrics/icv_values.csv'));

hcpd_master_table = innerjoin(hcpd_subs_master_table, hcpd_fd_table, 'Keys', 'src_subject_id');
hcpd_master_table = innerjoin(hcpd_master_table, hcpd_icv_table, 'Keys', 'src_subject_id');

hcpd_files = dir(fullfile(hcpd_input_dir, 'CorrVals_AmygHipOnly_*.mat'));
n_subj_hcpd = length(hcpd_files);
delta_r_hcpd = nan(n_subj_hcpd, n_vox_common);
subj_ids_hcpd = cell(n_subj_hcpd, 1);

for s = 1:n_subj_hcpd
    D = load(fullfile(hcpd_input_dir, hcpd_files(s).name));
    
    subj_ids_hcpd{s} = D.subj_id; 
    
    c_amyg = zeros(1, n_vox_common); c_hipp = zeros(1, n_vox_common);
    is_left = ismember(D.subcortical_labels, [9,10]); is_right = ismember(D.subcortical_labels, [1,2]);
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
sex_numeric(strcmpi(hcpd_final_covariates.sex, 'M')) = 1;
sex_numeric(strcmpi(hcpd_final_covariates.sex, 'F')) = 2;
hcpd_final_covariates.sex = sex_numeric;
[~,~,hcpd_final_covariates.race] = unique(hcpd_final_covariates.race);

cols_to_keep = {'src_subject_id', 'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV'};
master_covariates = hcpd_final_covariates(:, cols_to_keep);
master_covariates.interview_age = master_covariates.interview_age / 12;

[master_covariates, removed_rows_idx] = rmmissing(master_covariates);
master_delta_r = hcpd_final_delta_r;
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

output_raw = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_voxel_HCPD_only.csv');
writematrix(harmonized_dat', output_raw);

deltar_a_harmonized = mean(harmonized_delta_r_voxels(:, amyg_cols), 2, 'omitnan');
deltar_h_harmonized = mean(harmonized_delta_r_voxels(:, hipp_cols), 2, 'omitnan');

harmonized_master_table = master_covariates;
harmonized_master_table.deltar_a = deltar_a_harmonized;
harmonized_master_table.deltar_h = deltar_h_harmonized;

output_filename = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_table_HCPD_only.csv');
writetable(harmonized_master_table, output_filename);