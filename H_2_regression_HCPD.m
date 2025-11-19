clear; clc;

Top = '/data/project';
Path = '/data/project/HCPDfMRI/fmriresults01/';
addpath(genpath(fullfile(Top, 'tools')));

subs_table = readtable(fullfile(Top, 'HCPDfMRI/ndar_subject01.txt'));
subjects_to_process_info = subs_table;
SubDirs = dir(fullfile(Path, 'HCD*'));
subject_folders = {SubDirs.name};

subcortical_r = [1, 2]; subcortical_l = [9, 10];
medial = [42,57,59,73,74,75,76,77,78,79,80,81,85,86,88,104,106,109,180,181,182,195,196,222,237,239,253,254,255,256,257,258,259,260,261,265,266,268,284,286,289,360,361,362,375,376];
lateral = [26,27,28,60,84,87,93,94,96,97,98,100,101,102,103,105,107,112,113,186,187,206,207,208,240,264,267,273,274,277,278,280,281,282,283,285,287,292,293,366,367,82,83,89,90,91,92,95,99,108,110,114,262,263,269,270,271,272,275,276,279,288,290,294];
ROI_cortex = [medial, lateral];
cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
parc = cifti_label.cdata;
subcort_all_idx = find(parc >= 1 & parc <= 16);
amyg_hip_vox_idx = find(ismember(parc, [subcortical_l, subcortical_r]));
local_idx_into_subcort = ismember(subcort_all_idx, amyg_hip_vox_idx);
ctx_vertex_idx = find(ismember(parc, ROI_cortex)); 

HCPref = load(fullfile(Top, 'derivatives/FC_reference_HCP/FC_reference_HCP_scan1.mat'));
ref = HCPref.FC_mean; ref = ref(:, ctx_vertex_idx);
ref_parc_subcort = parc(subcort_all_idx);
ref_hip_l  = zscore(mean(ref(ismember(ref_parc_subcort, 9), :), 1));
ref_amyg_l = zscore(mean(ref(ismember(ref_parc_subcort, 10), :), 1));
ref_hip_r  = zscore(mean(ref(ismember(ref_parc_subcort, 1), :), 1));
ref_amyg_r = zscore(mean(ref(ismember(ref_parc_subcort, 2), :), 1));

subcort_labels_in_profile = parc(amyg_hip_vox_idx);
left_vox_idx = find(ismember(subcort_labels_in_profile, subcortical_l));
right_vox_idx = find(ismember(subcort_labels_in_profile, subcortical_r));

output_dir = fullfile(Top, 'derivatives/HCPD/tstat_maps_final');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
run_names = {'rfMRI_REST1_AP', 'rfMRI_REST1_PA', 'rfMRI_REST2_AP', 'rfMRI_REST2_PA', ...
             'rfMRI_REST1a_AP', 'rfMRI_REST1a_PA', 'rfMRI_REST1b_PA', ...
             'rfMRI_REST2a_AP', 'rfMRI_REST2a_PA', 'rfMRI_REST2b_AP'};

for s = 1:height(subjects_to_process_info)
    subj_id = subjects_to_process_info.src_subject_id{s};
    subj_folder_idx = contains(subject_folders, strrep(subj_id, 'HCD', ''));
    if ~any(subj_folder_idx), continue; end
    subj_folder = subject_folders{subj_folder_idx};
    
    all_cleaned_residuals = {};
    at_least_one_run_found = false;

    for r = 1:numel(run_names)
        run_name = run_names{r};
        base_path = fullfile(Path, subj_folder, 'MNINonLinear/Results/', run_name);
        ts_path = fullfile(base_path, [run_name '_Atlas_MSMAll_hp0_clean.dtseries.nii']);
        fd_path = fullfile(base_path, 'Movement_RelativeRMS.txt');
        motion_path = fullfile(base_path, 'Movement_Regressors.txt');

        if ~isfile(ts_path) || ~isfile(fd_path) || ~isfile(motion_path), continue; end
        at_least_one_run_found = true;
        
        fd_col = load(fd_path);
        if (sum(fd_col > 0.5) / numel(fd_col)) * 100 > 20, continue; end
        
        rest = cifti_read(ts_path);
        bold_data_full = rest.cdata';

        scrub_idx = find(fd_col > 0.5);
        frames_to_censor = [];
        for k = 1:length(scrub_idx), frames_to_censor = [frames_to_censor, (scrub_idx(k)-1):(scrub_idx(k)+2)]; end
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
        
        motion_12p = load(motion_path); motion_12p = motion_12p(5:end, :);
        motion_24p = [motion_12p, motion_12p.^2];
        
        X = [ones(size(motion_24p, 1), 1), motion_24p];
        beta = X \ bold_data;
        residuals = bold_data - X * beta;
        
        all_cleaned_residuals{end+1} = normalize(residuals);
    end
    
    if numel(all_cleaned_residuals) < 4, continue; end
    
    final_residuals_cat = cat(1, all_cleaned_residuals{:});
    A = final_residuals_cat';
    
    ts_sub = A(subcort_all_idx, :);
    ts_ctx = A(ctx_vertex_idx, :);
    FC_vox_profile_all = corr(ts_sub', ts_ctx');
    FC_vox_profile = FC_vox_profile_all(local_idx_into_subcort, :);
    
    is_bad_vertex = isnan(FC_vox_profile(1, :));
    FC_profile_clean = FC_vox_profile(:, ~is_bad_vertex);
    
    X_ref_L = [ones(sum(~is_bad_vertex), 1), ref_hip_l(~is_bad_vertex)', ref_amyg_l(~is_bad_vertex)'];
    X_ref_R = [ones(sum(~is_bad_vertex), 1), ref_hip_r(~is_bad_vertex)', ref_amyg_r(~is_bad_vertex)'];
    df_clean = sum(~is_bad_vertex) - 3;
    if df_clean <= 0, continue; end

    t_vals = zeros(size(FC_profile_clean, 1), 1);
    
    for i = 1:length(left_vox_idx)
        y = FC_profile_clean(left_vox_idx(i), :)';
        if all(y==0) || all(isnan(y)), continue; end
        y_z = zscore(y);
        [b, ~, ~, ~, stats] = regress(y_z, X_ref_L);
        beta_diff = b(3) - b(2);
        cov_b = inv(X_ref_L' * X_ref_L) * stats(4);
        se_diff = sqrt(cov_b(2,2) + cov_b(3,3) - 2 * cov_b(2,3));
        t_vals(left_vox_idx(i)) = beta_diff / se_diff;
    end

    for i = 1:length(right_vox_idx)
        y = FC_profile_clean(right_vox_idx(i), :)';
        if all(y==0) || all(isnan(y)), continue; end
        y_z = zscore(y);
        [b, ~, ~, ~, stats] = regress(y_z, X_ref_R);
        beta_diff = b(3) - b(2);
        cov_b = inv(X_ref_R' * X_ref_R) * stats(4);
        se_diff = sqrt(cov_b(2,2) + cov_b(3,3) - 2 * cov_b(2,3));
        t_vals(right_vox_idx(i)) = beta_diff / se_diff;
    end
    
    p_vals = 2 * tcdf(-abs(t_vals), df_clean);

    out_struct = struct();
    out_struct.subj_id = subj_id;
    out_struct.age = subjects_to_process_info.interview_age(s) / 12;
    out_struct.voxel_indices = amyg_hip_vox_idx;
    out_struct.subcortical_labels = parc(amyg_hip_vox_idx);
    out_struct.t_stat_map = t_vals;
    out_struct.p_vals = p_vals;
    
    save(fullfile(output_dir, sprintf('Tstat_AmygHipCompare_%s.mat', subj_id)), '-struct', 'out_struct');
end