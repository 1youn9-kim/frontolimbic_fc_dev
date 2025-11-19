clear; clc;


Top = '/data/project';
bcp_project_dir = fullfile(Top, 'BABY/image03/BCP');
fmriprep_dir = fullfile(bcp_project_dir, 'derivatives');
subject_list_file = fullfile(bcp_project_dir, 'final_sublist.txt');
addpath(genpath(fullfile(Top, 'tools')));
session_fc_output_dir = fullfile(Top, 'derivatives', 'BCP', 'session_level_hipp_fc_maps_final');
if ~exist(session_fc_output_dir, 'dir'), mkdir(session_fc_output_dir); end
hcpd_fc_output_dir = fullfile(Top, 'derivatives', 'HCPD', 'session_level_hipp_fc_maps');
if ~exist(hcpd_fc_output_dir, 'dir'), mkdir(hcpd_fc_output_dir); end
final_binned_output_dir = fullfile(Top, 'derivatives', 'BCP', 'binned_hipp_fc_maps_BCP');
if ~exist(final_binned_output_dir, 'dir'), mkdir(final_binned_output_dir); end
force_rerun_stage1 = false; 
force_rerun_hcpd_calc = false;

subcortical_hipp_r = 1; subcortical_hipp_l = 9;
medial = [42,57,59,73,74,75,76,77,78,79,80,81,85,86,88,104,106,109,180,181,182,195,196,222,237,239,253,254,255,256,257,258,259,260,261,265,266,268,284,286,289,360,361,362,375,376];
lateral = [26,27,28,60,84,87,93,94,96,97,98,100,101,102,103,105,107,112,113,186,187,206,207,208,240,264,267,273,274,277,278,280,281,282,283,285,287,292,293,366,367,82,83,89,90,91,92,95,99,108,110,114,262,263,269,270,271,272,275,276,279,288,290,294];
ROI_frontal_cortex = [medial, lateral];
atlas_path = fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
cifti_label = cifti_read(atlas_path); parc = cifti_label.cdata;
subcort_all_idx = find(parc >= 1 & parc <= 16);
hipp_mask_subcort = ismember(parc(subcort_all_idx), [subcortical_hipp_l, subcortical_hipp_r]);
ctx_vertex_idx = find(ismember(parc, ROI_frontal_cortex));
medial_idx = find(ismember(parc(ctx_vertex_idx), medial));
lateral_idx = find(ismember(parc(ctx_vertex_idx), lateral));

%% Calculate FC Maps
subjects_table = readtable(subject_list_file, 'ReadVariableNames', false);
subjects_to_process = cell(height(subjects_table), 1);
for i = 1:height(subjects_table)
    subjects_to_process{i} = sprintf('%06d', subjects_table.Var1(i));
end

session_files_exist = ~isempty(dir(fullfile(session_fc_output_dir, 'SessionFC_*.mat')));
if force_rerun_stage1 || ~session_files_exist
    dropout_log = {'SubjectID_Session', 'Reason'};
    fd_table_path = fullfile(Top, 'derivatives', 'BCP', 'mean_fd_run_BCP.csv');

    mean_fd_table = readtable(fd_table_path);

    for s = 1:numel(subjects_to_process)
        subj_id_short = subjects_to_process{s};
        subj_id = ['sub-' subj_id_short];
        subj_path = fullfile(fmriprep_dir, subj_id);
        if ~exist(subj_path, 'dir'), continue; end
        
        session_dirs = dir(fullfile(subj_path, 'ses-*'));
        for ses_idx = 1:numel(session_dirs)
            ses_name = session_dirs(ses_idx).name;
            
            sess_runs = mean_fd_table(strcmp(mean_fd_table.src_subject_id, subj_id) & strcmp(mean_fd_table.session, ses_name), :);
            if isempty(sess_runs)
                dropout_log(end+1,:) = {[subj_id '_' ses_name], 'No FD data found'};
                continue;
            end
            sess_runs = sess_runs(sess_runs.MeanFD <= 0.5, :);
            if isempty(sess_runs)
                 dropout_log(end+1,:) = {[subj_id '_' ses_name], 'All runs failed FD > 0.5'};
                 continue;
            end
            sess_runs = sortrows(sess_runs, 'MeanFD', 'ascend');
            best_run_idx = sess_runs.run(1); 
            run_str = num2str(best_run_idx);

            func_path = fullfile(subj_path, ses_name, 'func');
            best_run_file = dir(fullfile(func_path, [subj_id '_' ses_name '_*run-' run_str '*_space-fsLR_den-91k_bold.dtseries.nii']));
            if isempty(best_run_file)
                 run_str_pad = sprintf('%02d', best_run_idx);
                 best_run_file = dir(fullfile(func_path, [subj_id '_' ses_name '_*run-' run_str_pad '*_space-fsLR_den-91k_bold.dtseries.nii']));
            end
            if isempty(best_run_file)
                dropout_log(end+1,:) = {[subj_id '_' ses_name], sprintf('Best run (run-%s) file not found', run_str)}; 
                continue; 
            end
            ts_path = fullfile(best_run_file(1).folder, best_run_file(1).name);
            confounds_path = strrep(ts_path, '_space-fsLR_den-91k_bold.dtseries.nii', '_desc-confounds_timeseries.tsv');
            if ~isfile(ts_path) || ~isfile(confounds_path)
                 dropout_log(end+1,:) = {[subj_id '_' ses_name], 'Missing best run files'};
                 continue; 
            end

            confounds_table_full = readtable(confounds_path, 'FileType', 'text', 'Delimiter', '\t');
            fd_col = confounds_table_full.framewise_displacement;
            fd_col(isnan(fd_col)) = 0;
            if (sum(fd_col > 0.5) / numel(fd_col)) * 100 > 20
                dropout_log(end+1,:) = {[subj_id '_' ses_name], 'Best run failed QC'};
                continue;
            end
            bold_data_full = cifti_read(ts_path).cdata';
            scrub_idx = find(fd_col > 0.5);
            frames_to_censor = [];
            for k = 1:length(scrub_idx)
                frames_to_censor = [frames_to_censor, (scrub_idx(k)-1):(scrub_idx(k)+2)];
            end
            frames_to_censor = unique(frames_to_censor);
            frames_to_censor(frames_to_censor < 1 | frames_to_censor > size(bold_data_full, 1)) = [];
            bold_data_interp = bold_data_full;

            if ~isempty(frames_to_censor)
                good_frames = setdiff(1:size(bold_data_full, 1), frames_to_censor);
                for v = 1:size(bold_data_full, 2)
                    interpolated_vals = interp1(good_frames, bold_data_full(good_frames, v), frames_to_censor, 'linear', 'extrap');
                    bold_data_interp(frames_to_censor, v) = interpolated_vals;
                end
            end

            bold_data = bold_data_interp(6:end, :);
            confounds_table = confounds_table_full(6:end, :);
            mot6  = [confounds_table.trans_x, confounds_table.trans_y, confounds_table.trans_z, confounds_table.rot_x, confounds_table.rot_y, confounds_table.rot_z];
            dmot6 = [zeros(1,6); diff(mot6)];
            mot24 = [mot6, dmot6, mot6.^2, dmot6.^2];
            wm_csf_cols = {'white_matter', 'csf'};
            dct_cols = confounds_table.Properties.VariableNames(startsWith(confounds_table.Properties.VariableNames, 'cosine'));
            nuisance_matrix = [mot24, table2array(confounds_table(:, [wm_csf_cols, dct_cols]))];
            nuisance_matrix(isnan(nuisance_matrix)) = 0;
            X = [ones(size(nuisance_matrix, 1), 1), nuisance_matrix];
            beta = X \ bold_data;
            residuals = bold_data - X * beta;

            A = normalize(residuals', 2);
            ts_sub = A(subcort_all_idx, :); 
            ts_ctx = A(ctx_vertex_idx, :);

            FC_vox_profile = corr(ts_sub', ts_ctx');

            hipp_fc_patterns = FC_vox_profile(hipp_mask_subcort, :);
            session_mean_hipp_fc = mean(hipp_fc_patterns, 1, 'omitnan');
            
            out_struct = struct();
            out_struct.subj_id = subj_id_short;
            out_struct.ses_name = ses_name;
            out_struct.run_used = run_str;
            out_struct.mean_hipp_fc_vs_frontal_cortex = session_mean_hipp_fc;
            output_filename = fullfile(session_fc_output_dir, sprintf('SessionFC_sub-%s_%s.mat', subj_id_short, ses_name));
            save(output_filename, '-struct', 'out_struct');
        end
    end
    output_table = cell2table(dropout_log(2:end,:), 'VariableNames', dropout_log(1,:));
    writetable(output_table, fullfile(session_fc_output_dir, 'dropout_log.csv'));
else
    fprintf('files exist!\n');
end

%% Group Bin
main_demographics_file = fullfile(Top, 'BABY/image03/BCP/ndar_subject01.txt');
deltar_file = fullfile(Top, 'derivatives', 'MASTER_HARMONIZED_deltar_table_BCPHCPD.csv');
deltar_gamlss_file = fullfile(Top, 'derivatives', 'gamlss_scores_BCPHCPD', 'deltar_gamlss_zscores_BCPHCPD.csv');
bcp_icv_file = fullfile(Top, 'derivatives', 'BCP', 'icv_values_BCP.csv');
bcp_mean_fd_file = fullfile(Top, 'derivatives', 'BCP', 'mean_fd_run_BCP.csv');
main_demo_table = readtable(main_demographics_file);
deltar_table = readtable(deltar_file);
deltac_table = readtable(deltar_gamlss_file);
bcp_icv_table = readtable(bcp_icv_file);
bcp_mean_fd_table = readtable(bcp_mean_fd_file);

behavior_var = 'deltar_h_zscore';
main_demo_vars = {'src_subject_id', 'sex', 'race', 'site'};
deltar_vars = {'src_subject_id', 'interview_age'};
deltac_vars = {'src_subject_id', behavior_var};

subs_master_table = outerjoin(unique(main_demo_table(:, main_demo_vars),'rows'), unique(deltar_table(:, deltar_vars),'rows'), 'Keys', 'src_subject_id', 'MergeKeys', true);
subs_master_table = outerjoin(subs_master_table, unique(deltac_table(:, deltac_vars),'rows'), 'Keys', 'src_subject_id', 'MergeKeys', true);
subs_master_table = renamevars(subs_master_table, behavior_var, 'BehaviorZscore');

HCPref = load(fullfile(Top, 'derivatives/FC_reference_HCP/FC_reference_HCP_scan1.mat'));
ref = HCPref.FC_mean; 
ref = ref(:, ctx_vertex_idx);
ref_parc_subcort = parc(subcort_all_idx);

ref_hipp_l = mean(ref(ismember(ref_parc_subcort, 9), :), 1);
ref_hipp_r = mean(ref(ismember(ref_parc_subcort, 1), :), 1);

ref_vec = (ref_hipp_l + ref_hipp_r) / 2;

session_files = dir(fullfile(session_fc_output_dir, 'SessionFC_sub-*.mat'));
if isempty(session_files), error('No session FC files found.'); end
n_sessions = numel(session_files);
n_vertices = length(ctx_vertex_idx);
raw_fc_data = nan(n_sessions, n_vertices);
subject_ids_from_files = cell(n_sessions, 1);
session_names_from_files = cell(n_sessions, 1);
run_used_from_files = zeros(n_sessions, 1);
for i = 1:n_sessions
    data = load(fullfile(session_files(i).folder, session_files(i).name));
    if isfield(data, 'mean_hipp_fc_vs_frontal_cortex')
        raw_fc_data(i, :) = data.mean_hipp_fc_vs_frontal_cortex;
    elseif isfield(data, 'mean_amyg_fc_vs_frontal_cortex')
         raw_fc_data(i, :) = data.mean_amyg_fc_vs_frontal_cortex;
    end
    subject_ids_from_files{i} = ['sub-' data.subj_id];
    session_names_from_files{i} = data.ses_name;
    if isfield(data, 'run_used')
        run_used_from_files(i) = str2double(string(data.run_used));
    else
        run_used_from_files(i) = 1; 
    end
end

session_map = table(subject_ids_from_files, session_names_from_files, run_used_from_files, (1:n_sessions)', ...
    'VariableNames', {'src_subject_id', 'session', 'run', 'orig_row_idx'});

bcp_mean_fd_table.run = str2double(string(bcp_mean_fd_table.run));

session_map = innerjoin(session_map, bcp_mean_fd_table, 'Keys', {'src_subject_id', 'session', 'run'});
session_map = outerjoin(session_map, subs_master_table, 'Type', 'left', 'Keys', 'src_subject_id', 'MergeKeys', true);
session_map = outerjoin(session_map, unique(bcp_icv_table(:,{'src_subject_id','ICV'}),'rows'), 'Type', 'left', 'Keys', 'src_subject_id', 'MergeKeys', true);

required_vars = {'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV', 'BehaviorZscore'};
mask_complete = true(height(session_map), 1);
for v = 1:numel(required_vars)
    col = session_map.(required_vars{v});
    if isnumeric(col), mask_complete = mask_complete & ~isnan(col);
    else, mask_complete = mask_complete & ~ismissing(col); end
end
fprintf('Dropping %d sessions due to missing metadata.\n', sum(~mask_complete));
session_map = session_map(mask_complete, :);

session_map = sortrows(session_map, {'src_subject_id', 'interview_age'}, {'ascend', 'descend'});
[~, unique_idx] = unique(session_map.src_subject_id, 'first');
subs_metadata = session_map(unique_idx, :);

raw_fc_data = raw_fc_data(subs_metadata.orig_row_idx, :);

if ~iscell(subs_metadata.src_subject_id), subs_metadata.src_subject_id = cellstr(subs_metadata.src_subject_id); end


%% 5. Generate/Load Cached HCPD FC Maps
hcpd_ts_dir = fullfile(Top, 'HCPDfMRI/fmriresults01');
hcpd_subs_master_table = readtable(fullfile(Top, 'HCPDfMRI/ndar_subject01.txt'));
hcpd_fd_table = readtable(fullfile(Top, 'derivatives/qc_metrics/mean_fd_allsubs.csv'));
hcpd_icv_table = readtable(fullfile(Top, 'derivatives/qc_metrics/icv_values.csv'));
hcpd_gamlss_file = fullfile(Top, 'derivatives/gamlss_scores_multicohort_HCPD_BCP_last_final/deltac_gamlss_zscores_multicohort_HCPDBCP_last_CV_final_final.csv');
hcpd_gamlss_table = readtable(hcpd_gamlss_file);
hcpd_master_table = innerjoin(hcpd_subs_master_table, hcpd_fd_table, 'Keys', 'src_subject_id');
hcpd_master_table = innerjoin(hcpd_master_table, hcpd_icv_table, 'Keys', 'src_subject_id');
hcpd_master_table = innerjoin(hcpd_master_table, hcpd_gamlss_table, 'Keys', 'src_subject_id');
hcpd_master_table = renamevars(hcpd_master_table, 'deltar_h_zscore', 'BehaviorZscore');
all_fmri_files = dir(fullfile(hcpd_ts_dir, '*_MR'));
fmri_ids = extractBefore({all_fmri_files.name}, '_V1_MR')';
[~, subs_idx, ~] = intersect(hcpd_master_table.src_subject_id, fmri_ids);
hcpd_subs_to_process = hcpd_master_table(subs_idx, :);
n_subj_hcpd = height(hcpd_subs_to_process);
hipp_vox_idx_hcpd = find(ismember(parc, [1, 9])); 
for s = 1:n_subj_hcpd
    subj_id = hcpd_subs_to_process.src_subject_id{s};
    output_filename = fullfile(hcpd_fc_output_dir, sprintf('SessionFC_sub-%s.mat', subj_id));
    if force_rerun_hcpd_calc || ~exist(output_filename, 'file')
        ts_path = fullfile(hcpd_ts_dir, [subj_id '_V1_MR'], 'MNINonLinear/Results/rfMRI_REST/rfMRI_REST_Atlas_MSMAll_hp0_clean.dtseries.nii');
        if ~isfile(ts_path), continue; end
        rest = cifti_read(ts_path);
        A = normalize(rest.cdata, 2);
        ts_hipp_avg = mean(A(hipp_vox_idx_hcpd, :), 1);
        ts_ctx = A(ctx_vertex_idx, :);
        out_struct = struct();
        out_struct.subj_id = subj_id;
        out_struct.mean_hipp_fc_vs_frontal_cortex = corr(ts_hipp_avg', ts_ctx');
        save(output_filename, '-struct', 'out_struct');
    end
end
hcpd_files = dir(fullfile(hcpd_fc_output_dir, 'SessionFC_sub-*.mat'));
n_subj_hcpd = length(hcpd_files);
hcpd_raw_fc_data = nan(n_subj_hcpd, n_vertices);
subj_ids_hcpd = cell(n_subj_hcpd, 1);
for s = 1:n_subj_hcpd
    D = load(fullfile(hcpd_files(s).folder, hcpd_files(s).name));
    if isfield(D, 'mean_hipp_fc_vs_frontal_cortex')
        hcpd_raw_fc_data(s,:) = D.mean_hipp_fc_vs_frontal_cortex;
    elseif isfield(D, 'mean_amyg_fc_vs_frontal_cortex')
        hcpd_raw_fc_data(s,:) = D.mean_amyg_fc_vs_frontal_cortex;
    else
        error('File %s missing FC data field.', hcpd_files(s).name);
    end
    subj_ids_hcpd{s} = D.subj_id;
end
[~, master_idx] = ismember(subj_ids_hcpd, hcpd_master_table.src_subject_id);
aligned_cov_table = hcpd_master_table(master_idx(master_idx > 0), :);
aligned_fc_data = hcpd_raw_fc_data(master_idx > 0, :);
fd_mask = aligned_cov_table.MeanFD <= 0.2;
hcpd_final_covariates = aligned_cov_table(fd_mask, :);
hcpd_final_fc_data = aligned_fc_data(fd_mask, :);

%% 6. Harmonize with ComBat
bcp_final_covariates = subs_metadata;
bcp_final_fc_data = raw_fc_data;
sex_numeric_hcpd = nan(height(hcpd_final_covariates), 1);
sex_numeric_hcpd(strcmpi(hcpd_final_covariates.sex, 'M')) = 1;
sex_numeric_hcpd(strcmpi(hcpd_final_covariates.sex, 'F')) = 2;
hcpd_final_covariates.sex = sex_numeric_hcpd;
[~,~,hcpd_final_covariates.race] = unique(hcpd_final_covariates.race);
sex_numeric_bcp = nan(height(bcp_final_covariates), 1);
sex_numeric_bcp(strcmpi(bcp_final_covariates.sex, 'M')) = 1;
sex_numeric_bcp(strcmpi(bcp_final_covariates.sex, 'F')) = 2;
bcp_final_covariates.sex = sex_numeric_bcp;
[~,~,bcp_final_covariates.race] = unique(bcp_final_covariates.race);


cols_to_keep = {'src_subject_id', 'interview_age', 'sex', 'race', 'site', 'MeanFD', 'ICV', 'BehaviorZscore'};
hcpd_final_covariates.Dataset = repmat({'HCPD'}, height(hcpd_final_covariates), 1);
bcp_final_covariates.Dataset = repmat({'BCP'}, height(bcp_final_covariates), 1);
master_covariates = [hcpd_final_covariates(:, [cols_to_keep, 'Dataset']); bcp_final_covariates(:, [cols_to_keep, 'Dataset'])];
master_fc_data = [hcpd_final_fc_data; bcp_final_fc_data];
master_fc_data(master_fc_data > 1) = 1; master_fc_data(master_fc_data < -1) = -1;
[master_covariates, removed_rows_idx] = rmmissing(master_covariates);
master_fc_data(removed_rows_idx,:) = [];
master_fc_data_z = atanh(master_fc_data);
is_bcp_row = strcmp(master_covariates.Dataset, 'BCP');

dat_to_harmonize = master_fc_data_z';
dat_to_harmonize(isnan(dat_to_harmonize)) = 0;

[~, ~, batch] = unique(master_covariates.site);
age_z = zscore(master_covariates.interview_age);
sex_dummy = dummyvar(master_covariates.sex);
race_dummy = dummyvar(master_covariates.race);
meanfd_z = zscore(master_covariates.MeanFD);
icv_z = zscore(master_covariates.ICV);

mod = [age_z, sex_dummy(:, 1:end-1), race_dummy(:, 1:end-1), meanfd_z, icv_z];

harmonized_dat = combat(dat_to_harmonize, batch, mod, 1);
harmonized_master_fc_z = harmonized_dat';

raw_fc_data_z = harmonized_master_fc_z(is_bcp_row, :);
subs_metadata = master_covariates(is_bcp_row, :);

%% 7. Intra-session Normalization
mean_per_row = mean(raw_fc_data_z, 2, 'omitnan');
std_per_row  = std(raw_fc_data_z, 0, 2, 'omitnan');
std_per_row(std_per_row == 0) = 1; 
harmonized_data_z = (raw_fc_data_z - mean_per_row) ./ std_per_row;

%% Age Bins
age_bins = 0:1:5;
template_scalar = cifti_label;
template_scalar.cdata = nan(size(cifti_label.cdata,1), 1);
template_scalar.diminfo{2} = cifti_diminfo_make_scalars(1);
corrmat_medial = zeros(length(age_bins),3); corrmat_lateral = zeros(length(age_bins),3); corrmat_combined = zeros(length(age_bins),3);
corrmat_medial_all = zeros(length(age_bins),1); corrmat_lateral_all = zeros(length(age_bins),1);
for i = 1:length(age_bins)
    age_val = age_bins(i);
    age_mask = (subs_metadata.interview_age >= age_val) & (subs_metadata.interview_age < age_val + 1);
    if ~any(age_mask), continue; end
    sessions_in_age_bin_indices = find(age_mask);
    group_avg_all_z = mean(harmonized_data_z(sessions_in_age_bin_indices, :), 1, 'omitnan');
    final_group_map_all_r = tanh(group_avg_all_z);
    output_cifti_all = template_scalar;
    output_cifti_all.cdata(ctx_vertex_idx) = final_group_map_all_r';
    cifti_write(output_cifti_all, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_AgeBin%02d_all.dscalar.nii', age_val)));
    corrmat_medial_all(i) = corr(ref_vec(medial_idx)', final_group_map_all_r(medial_idx)', 'rows','complete');
    corrmat_lateral_all(i) = corr(ref_vec(lateral_idx)', final_group_map_all_r(lateral_idx)', 'rows','complete');
    z_scores_in_bin = subs_metadata.BehaviorZscore(sessions_in_age_bin_indices);
    if numel(z_scores_in_bin) < 3, continue; end
    behavior_tertiles = quantile(z_scores_in_bin, [1/3, 2/3]);
    group_indices = struct('low', sessions_in_age_bin_indices(z_scores_in_bin < behavior_tertiles(1)), ...
                           'mid', sessions_in_age_bin_indices(z_scores_in_bin >= behavior_tertiles(1) & z_scores_in_bin < behavior_tertiles(2)), ...
                           'high', sessions_in_age_bin_indices(z_scores_in_bin >= behavior_tertiles(2)));
    group_keys = {'low', 'mid', 'high'}; group_labels = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    for g = 1:length(group_keys)
        idx = group_indices.(group_keys{g});
        if isempty(idx), continue; end
        group_avg_z = mean(harmonized_data_z(idx, :), 1, 'omitnan');
        final_group_map_r = tanh(group_avg_z);
        output_cifti = template_scalar;
        output_cifti.cdata(ctx_vertex_idx) = final_group_map_r';
        cifti_write(output_cifti, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_AgeBin%02d_%s.dscalar.nii', age_val, group_labels{g})));
        corrmat_medial(i,g) = corr(ref_vec(medial_idx)', final_group_map_r(medial_idx)', 'rows','complete');
        corrmat_lateral(i,g) = corr(ref_vec(lateral_idx)', final_group_map_r(lateral_idx)', 'rows','complete');
        corrmat_combined(i,g) = corr(ref_vec', final_group_map_r', 'rows','complete');
    end
end
writematrix(corrmat_medial, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBinYearly_CorrWRef_medial.csv'));
writematrix(corrmat_lateral, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBinYearly_CorrWRef_lateral.csv'));
writematrix(corrmat_combined, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBinYearly_CorrWRef.csv'));
writematrix(corrmat_medial_all, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBinYearly_CorrWRef_medial_all.csv'));
writematrix(corrmat_lateral_all, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBinYearly_CorrWRef_lateral_all.csv'));

%% Tertiles
subject_ages = subs_metadata.interview_age;
age_tertiles = quantile(subject_ages, [1/3, 2/3]);
age_bin_masks = {subject_ages < age_tertiles(1), subject_ages >= age_tertiles(1) & subject_ages < age_tertiles(2), subject_ages >= age_tertiles(2)};
age_bin_labels = {'Age_Tertile1', 'Age_Tertile2', 'Age_Tertile3'};
corrmat_medial_3split = zeros(3,3); corrmat_lateral_3split = zeros(3,3); corrmat_combined_3split = zeros(3,3);
corrmat_medial_all_3split = zeros(3,1); corrmat_lateral_all_3split = zeros(3,1);
for i = 1:3
    mask = age_bin_masks{i};
    if ~any(mask), continue; end
    idx = find(mask);
    group_avg_all_z = mean(harmonized_data_z(idx, :), 1, 'omitnan');
    final_group_map_all_r = tanh(group_avg_all_z);
    output_cifti_all = template_scalar;
    output_cifti_all.cdata(ctx_vertex_idx) = final_group_map_all_r';
    cifti_write(output_cifti_all, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_%s_all.dscalar.nii', age_bin_labels{i})));
    corrmat_medial_all_3split(i) = corr(ref_vec(medial_idx)', final_group_map_all_r(medial_idx)', 'rows','complete');
    corrmat_lateral_all_3split(i) = corr(ref_vec(lateral_idx)', final_group_map_all_r(lateral_idx)', 'rows','complete');
    z_scores = subs_metadata.BehaviorZscore(idx);
    beh_tertiles = quantile(z_scores, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scores < beh_tertiles(1)), 'mid', idx(z_scores >= beh_tertiles(1) & z_scores < beh_tertiles(2)), 'high', idx(z_scores >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labels = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    for g = 1:length(g_keys)
        s_idx = g_idx.(g_keys{g});
        if isempty(s_idx), continue; end
        group_avg_z = mean(harmonized_data_z(s_idx, :), 1, 'omitnan');
        final_group_map_r = tanh(group_avg_z);
        output_cifti = template_scalar;
        output_cifti.cdata(ctx_vertex_idx) = final_group_map_r';
        cifti_write(output_cifti, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_%s_%s.dscalar.nii', age_bin_labels{i}, g_labels{g})));
        corrmat_medial_3split(i,g) = corr(ref_vec(medial_idx)', final_group_map_r(medial_idx)', 'rows','complete');
        corrmat_lateral_3split(i,g) = corr(ref_vec(lateral_idx)', final_group_map_r(lateral_idx)', 'rows','complete');
        corrmat_combined_3split(i,g) = corr(ref_vec', final_group_map_r', 'rows','complete');
    end
end
writematrix(corrmat_medial_3split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin3Split_CorrWRef_medial.csv'));
writematrix(corrmat_lateral_3split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin3Split_CorrWRef_lateral.csv'));
writematrix(corrmat_combined_3split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin3Split_CorrWRef.csv'));
writematrix(corrmat_medial_all_3split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin3Split_CorrWRef_medial_all.csv'));
writematrix(corrmat_lateral_all_3split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin3Split_CorrWRef_lateral_all.csv'));

%% Quartiles
age_quartiles = quantile(subject_ages, [0.25, 0.5, 0.75]);
age_bin_masks = {subject_ages < age_quartiles(1), subject_ages >= age_quartiles(1) & subject_ages < age_quartiles(2), subject_ages >= age_quartiles(2) & subject_ages < age_quartiles(3), subject_ages >= age_quartiles(3)};
age_bin_labels = {'Age_Quartile1', 'Age_Quartile2', 'Age_Quartile3', 'Age_Quartile4'};
corrmat_medial_4split = zeros(4,3); corrmat_lateral_4split = zeros(4,3); corrmat_combined_4split = zeros(4,3);
corrmat_medial_all_4split = zeros(4,1); corrmat_lateral_all_4split = zeros(4,1);
for i = 1:4
    mask = age_bin_masks{i};
    if ~any(mask), continue; end
    idx = find(mask);
    group_avg_all_z = mean(harmonized_data_z(idx, :), 1, 'omitnan');
    final_group_map_all_r = tanh(group_avg_all_z);
    output_cifti_all = template_scalar;
    output_cifti_all.cdata(ctx_vertex_idx) = final_group_map_all_r';
    cifti_write(output_cifti_all, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_%s_all.dscalar.nii', age_bin_labels{i})));
    corrmat_medial_all_4split(i) = corr(ref_vec(medial_idx)', final_group_map_all_r(medial_idx)', 'rows','complete');
    corrmat_lateral_all_4split(i) = corr(ref_vec(lateral_idx)', final_group_map_all_r(lateral_idx)', 'rows','complete');
    z_scores = subs_metadata.BehaviorZscore(idx);
    beh_tertiles = quantile(z_scores, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scores < beh_tertiles(1)), 'mid', idx(z_scores >= beh_tertiles(1) & z_scores < beh_tertiles(2)), 'high', idx(z_scores >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labels = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    for g = 1:length(g_keys)
        s_idx = g_idx.(g_keys{g});
        if isempty(s_idx), continue; end
        group_avg_z = mean(harmonized_data_z(s_idx, :), 1, 'omitnan');
        final_group_map_r = tanh(group_avg_z);
        output_cifti = template_scalar;
        output_cifti.cdata(ctx_vertex_idx) = final_group_map_r';
        cifti_write(output_cifti, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_%s_%s.dscalar.nii', age_bin_labels{i}, g_labels{g})));
        corrmat_medial_4split(i,g) = corr(ref_vec(medial_idx)', final_group_map_r(medial_idx)', 'rows','complete');
        corrmat_lateral_4split(i,g) = corr(ref_vec(lateral_idx)', final_group_map_r(lateral_idx)', 'rows','complete');
        corrmat_combined_4split(i,g) = corr(ref_vec', final_group_map_r', 'rows','complete');
    end
end
writematrix(corrmat_medial_4split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin4Split_CorrWRef_medial.csv'));
writematrix(corrmat_lateral_4split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin4Split_CorrWRef_lateral.csv'));
writematrix(corrmat_combined_4split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin4Split_CorrWRef.csv'));
writematrix(corrmat_medial_all_4split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin4Split_CorrWRef_medial_all.csv'));
writematrix(corrmat_lateral_all_4split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin4Split_CorrWRef_lateral_all.csv'));

%% Median Split
age_median = quantile(subject_ages, 0.5);
age_bin_masks = {subject_ages < age_median, subject_ages >= age_median};
age_bin_labels = {'Age_MedianSplit_Low', 'Age_MedianSplit_High'};
corrmat_medial_2split = zeros(2,3); corrmat_lateral_2split = zeros(2,3); corrmat_combined_2split = zeros(2,3);
corrmat_medial_all_2split = zeros(2,1); corrmat_lateral_all_2split = zeros(2,1);
for i = 1:2
    mask = age_bin_masks{i};
    if ~any(mask), continue; end
    idx = find(mask);
    group_avg_all_z = mean(harmonized_data_z(idx, :), 1, 'omitnan');
    final_group_map_all_r = tanh(group_avg_all_z);
    output_cifti_all = template_scalar;
    output_cifti_all.cdata(ctx_vertex_idx) = final_group_map_all_r';
    cifti_write(output_cifti_all, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_%s_all.dscalar.nii', age_bin_labels{i})));
    corrmat_medial_all_2split(i) = corr(ref_vec(medial_idx)', final_group_map_all_r(medial_idx)', 'rows','complete');
    corrmat_lateral_all_2split(i) = corr(ref_vec(lateral_idx)', final_group_map_all_r(lateral_idx)', 'rows','complete');
    z_scores = subs_metadata.BehaviorZscore(idx);
    beh_tertiles = quantile(z_scores, [1/3, 2/3]);
    g_idx = struct('low', idx(z_scores < beh_tertiles(1)), 'mid', idx(z_scores >= beh_tertiles(1) & z_scores < beh_tertiles(2)), 'high', idx(z_scores >= beh_tertiles(2)));
    g_keys = {'low', 'mid', 'high'}; g_labels = {'Bottom33pct', 'Middle33pct', 'Top33pct'};
    for g = 1:length(g_keys)
        s_idx = g_idx.(g_keys{g});
        if isempty(s_idx), continue; end
        group_avg_z = mean(harmonized_data_z(s_idx, :), 1, 'omitnan');
        final_group_map_r = tanh(group_avg_z);
        output_cifti = template_scalar;
        output_cifti.cdata(ctx_vertex_idx) = final_group_map_r';
        cifti_write(output_cifti, fullfile(final_binned_output_dir, sprintf('HippAvgFC_BCP_Denoised_%s_%s.dscalar.nii', age_bin_labels{i}, g_labels{g})));
        corrmat_medial_2split(i,g) = corr(ref_vec(medial_idx)', final_group_map_r(medial_idx)', 'rows','complete');
        corrmat_lateral_2split(i,g) = corr(ref_vec(lateral_idx)', final_group_map_r(lateral_idx)', 'rows','complete');
        corrmat_combined_2split(i,g) = corr(ref_vec', final_group_map_r', 'rows','complete');
    end
end
writematrix(corrmat_medial_2split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin2Split_CorrWRef_medial.csv'));
writematrix(corrmat_lateral_2split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin2Split_CorrWRef_lateral.csv'));
writematrix(corrmat_combined_2split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin2Split_CorrWRef.csv'));
writematrix(corrmat_medial_all_2split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin2Split_CorrWRef_medial_all.csv'));
writematrix(corrmat_lateral_all_2split, fullfile(final_binned_output_dir, 'HippAvgFC_BCP_Denoised_AgeBin2Split_CorrWRef_lateral_all.csv'));
