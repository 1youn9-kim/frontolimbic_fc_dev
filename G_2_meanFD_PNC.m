clear; clc; 

Top = '/data/project';
addpath(genpath(fullfile(Top, 'tools')));

input_dir = '/data/project/derivatives/PNC/corr_maps';
files = dir(fullfile(input_dir, 'CorrVals_AmygHipOnly_*.mat'));
tmp = load(fullfile(input_dir, files(1).name));
common_voxel_labels = tmp.subcortical_labels;
amyg_cols = ismember(common_voxel_labels, [2, 10]); 
hipp_cols = ismember(common_voxel_labels, [1, 9]);  
n_vox_common = length(common_voxel_labels);

subs_master_table = readtable(fullfile(Top, 'PNC_subjectlist.csv'));
subs_master_table.src_subject_id = [];
subs_master_table = renamevars(subs_master_table,'bids','src_subject_id');
fd_table = readtable(fullfile(Top, 'derivatives/qc_metrics/mean_fd_PNC.csv'));
icv_table = readtable(fullfile(Top, 'derivatives/qc_metrics/icv_values_PNC.csv'));
master_table = innerjoin(subs_master_table, fd_table, 'Keys', 'src_subject_id');
master_table = innerjoin(master_table, icv_table, 'Keys', 'src_subject_id');

n_subj = length(files);
delta_r = nan(n_subj, n_vox_common);
subj_ids = cell(n_subj, 1);
for s = 1:n_subj
    D = load(fullfile(input_dir, files(s).name));
    if ~isequal(D.subcortical_labels, common_voxel_labels), error('Voxel space mismatch!'); end
    subj_ids{s} = D.subj_id; 

    c_amyg = zeros(1, n_vox_common);
    c_hipp = zeros(1, n_vox_common);
    is_left = ismember(D.subcortical_labels, [9,10]);
    is_right = ismember(D.subcortical_labels, [1,2]);
    c_amyg(is_left)  = D.corr_left(is_left,2);
    c_hipp(is_left)  = D.corr_left(is_left,1);
    c_amyg(is_right) = D.corr_right(is_right,2);
    c_hipp(is_right) = D.corr_right(is_right,1);
    delta_r(s,:) = c_amyg - c_hipp;
end

[is_present, location_in_loaded_data] = ismember(master_table.src_subject_id, subj_ids);
aligned_cov_table = master_table(is_present, :);
reorder_indices = location_in_loaded_data(is_present);
aligned_delta_r = delta_r(reorder_indices, :);

pnc_final_covariates = aligned_cov_table;
pnc_final_delta_r = aligned_delta_r;

sex_numeric = nan(height(pnc_final_covariates), 1);
sex_numeric(strcmpi(pnc_final_covariates.sex, 'M')) = 1;
sex_numeric(strcmpi(pnc_final_covariates.sex, 'F')) = 2;
pnc_final_covariates.sex = sex_numeric;
[~,~,pnc_final_covariates.race] = unique(pnc_final_covariates.race);

cols_to_keep = {'src_subject_id', 'age', 'sex', 'race', 'site', 'MeanFD', 'ICV'};
master_covariates = pnc_final_covariates(:, cols_to_keep);
master_covariates.age = master_covariates.age / 12;

[master_covariates, removed_rows_idx] = rmmissing(master_covariates);
master_delta_r = pnc_final_delta_r;
master_delta_r(removed_rows_idx,:) = [];

deltar_a = mean(master_delta_r(:, amyg_cols), 2, 'omitnan');
deltar_h = mean(master_delta_r(:, hipp_cols), 2, 'omitnan');

pnc_master_table = master_covariates;
pnc_master_table.deltar_a = deltar_a;
pnc_master_table.deltar_h = deltar_h;

output_filename = fullfile(Top, 'derivatives', 'MASTER_deltar_table_PNC.csv');
writetable(pnc_master_table, output_filename);