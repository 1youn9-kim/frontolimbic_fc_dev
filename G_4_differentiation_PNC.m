
clear; clc;


Top = '/data/projects';
input_dir = '/data/projects/punim2400/derivatives/PNC/corr_maps';
out_dir = fullfile(input_dir, 'differentiation');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
addpath(genpath(fullfile(Top, 'tools')));

files = dir(fullfile(input_dir, 'CorrVals_AmygHipOnly_*.mat'));
n_subj = length(files);

tmp = load(fullfile(input_dir, files(1).name));
voxel_idx = tmp.voxel_indices;
n_vox = length(voxel_idx);

subs_master_table = readtable(fullfile(Top, 'PNC_subjectlist.csv'));
subs_master_table.src_subject_id = [];
subs_master_table = renamevars(subs_master_table,'bids','src_subject_id');
fd_table = readtable(fullfile(Top, 'derivatives/qc_metrics/mean_fd_PNC.csv'));
icv_table = readtable(fullfile(Top, 'derivatives/qc_metrics/icv_values_PNC.csv'));
master_table = innerjoin(subs_master_table, fd_table, 'Keys', 'src_subject_id');
master_table = innerjoin(master_table, icv_table, 'Keys', 'src_subject_id');

delta_r = nan(n_subj, n_vox);
subj_ids = cell(n_subj, 1);
for s = 1:n_subj
    D = load(fullfile(input_dir, files(s).name));
    subj_ids{s} = D.subj_id;

    c_amyg = zeros(1, n_vox);
    c_hipp = zeros(1, n_vox);
    is_left = ismember(D.subcortical_labels, [9,10]);
    is_right = ismember(D.subcortical_labels, [1,2]);
    c_amyg(is_left)  = D.corr_left(is_left,2);
    c_hipp(is_left)  = D.corr_left(is_left,1);
    c_amyg(is_right) = D.corr_right(is_right,2);
    c_hipp(is_right) = D.corr_right(is_right,1);
    delta_r(s,:) = atanh(c_amyg) - atanh(c_hipp);
end

is_missing_from_load = all(isnan(delta_r), 2);
loaded_subj_ids = subj_ids(~is_missing_from_load);
loaded_delta_r = delta_r(~is_missing_from_load, :);

[is_present, location_in_loaded_data] = ismember(master_table.src_subject_id, loaded_subj_ids);

final_cov_table = master_table(is_present, :);

reorder_indices = location_in_loaded_data(is_present);

final_delta_r = loaded_delta_r(reorder_indices, :);
final_subj_ids = loaded_subj_ids(reorder_indices, :);

save(fullfile(out_dir,'deltar_PNC.mat'),"final_delta_r");

age = final_cov_table.age / 12;
sex = categorical(final_cov_table.sex);
mean_fd = final_cov_table.MeanFD;
icv = final_cov_table.ICV;
race = categorical(final_cov_table.race);

X = [age, dummyvar(sex), dummyvar(race), mean_fd, icv];
X = normalize(X);

beta_age = nan(n_vox,1);
pvals_age = nan(n_vox,1);
for v = 1:n_vox
    y = final_delta_r(:,v);
    mdl = fitlm(X, y);
    beta_age(v) = mdl.Coefficients.Estimate(2);
    pvals_age(v) = mdl.Coefficients.pValue(2);
end

[~, ~, ~, qvals_age] = fdr_bh(pvals_age);

save(fullfile(out_dir, 'deltar_voxelwise_regression_PNC.mat'), ...
    'voxel_idx', 'final_subj_ids', 'beta_age', 'pvals_age', 'qvals_age');

template = cifti_read('/data/projects/punim2400/Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
template_scalar = template;
template_scalar.cdata = nan(size(template.cdata,1), 1);
template_scalar.diminfo{2} = cifti_diminfo_make_scalars(1);

% Beta map
beta_map = template_scalar;
beta_map.cdata(voxel_idx) = beta_age;
cifti_write(beta_map, fullfile(out_dir, 'deltac_beta1_standardized_robust_PNC.dscalar.nii'));

% Raw p
pval_map = template_scalar;
pval_map.cdata(voxel_idx) = pvals_age;
cifti_write(pval_map, fullfile(out_dir, 'deltac_pvals_standardized_robust_PNC.dscalar.nii'));

% FDR q
qval_map = template_scalar;
qval_map.cdata(voxel_idx) = qvals_age;
cifti_write(qval_map, fullfile(out_dir, 'deltac_qvals_standardized_robust_PNC.dscalar.nii'));