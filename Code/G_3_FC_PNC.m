clear; clc;


Top = '/data/project';
Path = fullfile(Top, 'PNC/Data');
addpath(genpath(fullfile(Top, 'tools')));

subcortical_r = [1, 2]; subcortical_l = [9, 10];
medial = [42,57,59,73,74,75,76,77,78,79,80,81,85,86,88,104,106,109,180,181,182,195,196,222,237,239,253,254,255,256,257,258,259,260,261,265,266,268,284,286,289,360,361,362,375,376];
lateral = [26,27,28,60,84,87,93,94,96,97,98,100,101,102,103,105,107,112,113,186,187,206,207,208,240,264,267,273,274,277,278,280,281,282,283,285,287,292,293,366,367,82,83,89,90,91,92,95,99,108,110,114,262,263,269,270,271,272,275,276,279,288,290,294];
ROI_cortex = [medial, lateral];
atlas_path = fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
cifti_label = cifti_read(atlas_path); parc = cifti_label.cdata;
subcort_all_idx = find(parc >= 1 & parc <= 16);
amyg_hip_vox_idx = find(ismember(parc, [subcortical_l, subcortical_r]));
local_idx_into_subcort = ismember(subcort_all_idx, amyg_hip_vox_idx);
ctx_vertex_idx = find(ismember(parc, ROI_cortex));

keeplist_file = fullfile(Top, 'PNC_subjectlist.csv');
keeplist_table = readtable(keeplist_file, 'ReadVariableNames', false);
subjects_to_keep_from_list = keeplist_table.Var1;
if isnumeric(subjects_to_keep_from_list), subjects_to_keep_from_list = cellstr(num2str(subjects_to_keep_from_list)); end
listing = dir(Path);
all_folders_on_disk = {listing([listing.isdir]).name};
subjects_on_disk = setdiff(all_folders_on_disk, {'.', '..'});
initial_subjects = intersect(subjects_on_disk, subjects_to_keep_from_list);

fd_threshold_qc = 0.5;
percent_threshold = 20;
is_subject_good = true(numel(initial_subjects), 1);

for s = 1:numel(initial_subjects)
    subj_id = initial_subjects{s};
    confounds_filename = sprintf('%s_task-rest_desc-confounds_timeseries.tsv', subj_id);
    confounds_path = fullfile(Path, subj_id, 'func', confounds_filename);
    if ~isfile(confounds_path)
        is_subject_good(s) = false;
        continue;
    end
    T = readtable(confounds_path, 'FileType', 'text', 'Delimiter', '\t');
    fd_col = T.framewise_displacement;
    percent_high_motion = (sum(fd_col > fd_threshold_qc, 'omitnan') / numel(fd_col)) * 100;
    if percent_high_motion > percent_threshold
        is_subject_good(s) = false;
    end
end
subjects_to_process = initial_subjects(is_subject_good);

HCPref = load(fullfile(Top, 'derivatives/FC_reference_HCP/FC_reference_HCP_scan1.mat'));
ref = HCPref.FC_mean; ref_ctx = ref(:, ctx_vertex_idx);
ref_parc_subcort = parc(subcort_all_idx);
ref_hip_l  = mean(ref_ctx(ismember(ref_parc_subcort, 9), :), 1);
ref_amyg_l = mean(ref_ctx(ismember(ref_parc_subcort, 10), :), 1);
ref_hip_r  = mean(ref_ctx(ismember(ref_parc_subcort, 1), :), 1);
ref_amyg_r = mean(ref_ctx(ismember(ref_parc_subcort, 2), :), 1);

output_dir = fullfile(Top, 'derivatives/PNC/corr_maps');

for s = 1:numel(subjects_to_process)
    subj_id = subjects_to_process{s};
    ts_filename = sprintf('%s_task-rest_space-fsLR_den-91k_bold.dtseries.nii', subj_id);
    ts_path = fullfile(Path, subj_id, 'func', ts_filename);
    if ~isfile(ts_path), continue; end

    rest = cifti_read(ts_path);
    bold_data_full = rest.cdata';
    confounds_filename = sprintf('%s_task-rest_desc-confounds_timeseries.tsv', subj_id);
    confounds_path = fullfile(Path, subj_id, 'func', confounds_filename);
    if ~isfile(confounds_path), continue; end
    confounds_table_full = readtable(confounds_path, 'FileType', 'text', 'Delimiter', '\t');
    
    fd_vals_full = confounds_table_full.framewise_displacement;
    fd_threshold_scrub = 0.5;
    scrub_idx = find(fd_vals_full > fd_threshold_scrub);
    
    frames_to_censor = [];
    for k = 1:length(scrub_idx)
        censor_win = (scrub_idx(k)-1):(scrub_idx(k)+2);
        frames_to_censor = [frames_to_censor, censor_win];
    end
    frames_to_censor = unique(frames_to_censor);
    frames_to_censor(frames_to_censor < 1 | frames_to_censor > size(bold_data_full, 1)) = [];

    bold_data_interp = bold_data_full;
    if ~isempty(frames_to_censor)
        good_frames = setdiff(1:size(bold_data_full, 1), frames_to_censor);
        for v = 1:size(bold_data_full, 2)
            vertex_ts = bold_data_full(:, v);
            interpolated_vals = interp1(good_frames, vertex_ts(good_frames), frames_to_censor, 'linear', 'extrap');
            bold_data_interp(frames_to_censor, v) = interpolated_vals;
        end
    end
    
    bold_data = bold_data_interp(5:end, :);
    confounds_table = confounds_table_full(5:end, :);

    motion_base = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
    motion_deriv = cellfun(@(x) [x, '_derivative1'], motion_base, 'UniformOutput', false);
    motion_24p = [motion_base, motion_deriv, cellfun(@(x) [x, '_power2'], motion_base, 'UniformOutput', false), cellfun(@(x) [x, '_power2'], motion_deriv, 'UniformOutput', false)];
    wm_csf = {'white_matter', 'csf'};
    dct_cols = confounds_table.Properties.VariableNames(startsWith(confounds_table.Properties.VariableNames, 'cosine'));
    
    nuisance_cols = [motion_24p, wm_csf, dct_cols];
    nuisance_matrix = table2array(confounds_table(:, nuisance_cols));
    nuisance_matrix(isnan(nuisance_matrix)) = 0;
    
    X = [ones(size(nuisance_matrix, 1), 1), nuisance_matrix];
    beta = X \ bold_data;
    residuals = bold_data - X * beta;
    
    A = normalize(residuals', 2); 
    
    ts_sub = A(subcort_all_idx, :); ts_ctx = A(ctx_vertex_idx, :);
    FC_vox_profile_all = corr(ts_sub', ts_ctx');
    FC_vox_profile = FC_vox_profile_all(local_idx_into_subcort, :);
    FC_vox_profile = atanh(FC_vox_profile);
    n_vox = size(FC_vox_profile, 1);
    corr_vals_left = zeros(n_vox, 2); corr_vals_right = zeros(n_vox, 2);
    for i = 1:n_vox
        corr_vals_left(i,1)  = corr(FC_vox_profile(i,:)', ref_hip_l', 'Rows', 'complete');
        corr_vals_left(i,2)  = corr(FC_vox_profile(i,:)', ref_amyg_l', 'Rows', 'complete');
        corr_vals_right(i,1) = corr(FC_vox_profile(i,:)', ref_hip_r', 'Rows', 'complete');
        corr_vals_right(i,2) = corr(FC_vox_profile(i,:)', ref_amyg_r', 'Rows', 'complete');
    end
    out_struct = struct();
    out_struct.subj_id = subj_id; out_struct.voxel_indices = amyg_hip_vox_idx;
    out_struct.subcortical_labels = parc(amyg_hip_vox_idx);
    out_struct.corr_left = corr_vals_left; out_struct.corr_right = corr_vals_right;
    output_filename = fullfile(output_dir, sprintf('CorrVals_AmygHipOnly_%s.mat', subj_id));
    save(output_filename, '-struct', 'out_struct');
end
