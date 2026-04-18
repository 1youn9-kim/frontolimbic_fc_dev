clear; clc;
Top = '/data/projects';
bcp_project_dir = fullfile(Top, 'BABY/image03/BCP');
nibabies_dir = fullfile(bcp_project_dir, 'derivatives', 'nibabies');
subject_list_file = fullfile(bcp_project_dir, 'final_sublist.txt');
output_dir = fullfile(Top, 'derivatives', 'BCP', 'corr_maps');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
addpath(genpath(fullfile(bcp_project_dir, 'tools')));

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

HCPref = load(fullfile(Top, 'derivatives/FC_reference_HCP/FC_reference_HCP_scan1.mat'));
ref = HCPref.FC_mean; ref_ctx = ref(:, ctx_vertex_idx);
ref_parc_subcort = parc(subcort_all_idx);
ref_hip_l  = mean(ref_ctx(ismember(ref_parc_subcort, 9), :), 1);
ref_amyg_l = mean(ref_ctx(ismember(ref_parc_subcort, 10), :), 1);
ref_hip_r  = mean(ref_ctx(ismember(ref_parc_subcort, 1), :), 1);
ref_amyg_r = mean(ref_ctx(ismember(ref_parc_subcort, 2), :), 1);

subjects_table = readtable(subject_list_file, 'ReadVariableNames', false);
subjects_to_process = cell(height(subjects_table), 1);
for i = 1:height(subjects_table), subjects_to_process{i} = sprintf('%06d', subjects_table.Var1(i)); end

dropout_log = {'SubjectID_Session', 'Reason'};
for s = 1:numel(subjects_to_process)
    subj_id = subjects_to_process{s};
    subj_path = fullfile(nibabies_dir, ['sub-' subj_id]);
    if ~exist(subj_path, 'dir'), continue; end
    
    session_dirs = dir(fullfile(subj_path, 'ses-*'));
    for ses_idx = 1:numel(session_dirs)
        ses_name = session_dirs(ses_idx).name;
        func_path = fullfile(subj_path, ses_name, 'func');
        if ~exist(func_path, 'dir'), continue; end
        bold_files = dir(fullfile(func_path, ['sub-' subj_id '_' ses_name '_*_space-fsLR_den-91k_bold.dtseries.nii']));
        if isempty(bold_files), dropout_log(end+1,:) = {[subj_id '_' ses_name], 'No dtseries files found'}; continue; end
        
        for r = 1:numel(bold_files)
            run = string(r); % Use sequential counter to match CSV
            ts_path = fullfile(bold_files(r).folder, bold_files(r).name);
            confounds_path = strrep(ts_path, '_space-fsLR_den-91k_bold.dtseries.nii', '_desc-confounds_timeseries.tsv');
            if ~isfile(ts_path) || ~isfile(confounds_path), continue; end
            
            confounds_table_full = readtable(confounds_path, 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'n/a','NA'});
            fd_col = confounds_table_full.framewise_displacement; fd_col(isnan(fd_col)) = 0;
            if (sum(fd_col > 0.5) / numel(fd_col)) * 100 > 20, continue; end
            
            bold_data_full = cifti_read(ts_path).cdata';
            scrub_idx = find(fd_col > 0.5);
            frames_to_censor = unique([scrub_idx-1; scrub_idx; scrub_idx+1; scrub_idx+2]);
            frames_to_censor(frames_to_censor < 1 | frames_to_censor > size(bold_data_full, 1)) = [];
            
            if ~isempty(frames_to_censor)
                good_frames = setdiff(1:size(bold_data_full, 1), frames_to_censor);
                for v = 1:size(bold_data_full, 2)
                    bold_data_full(frames_to_censor, v) = interp1(good_frames, bold_data_full(good_frames, v), frames_to_censor, 'linear', 'extrap');
                end
            end
            
            bold_data = bold_data_full(5:end, :);
            confounds_table = confounds_table_full(5:end, :);
            mot6  = [confounds_table.trans_x, confounds_table.trans_y, confounds_table.trans_z, confounds_table.rot_x, confounds_table.rot_y, confounds_table.rot_z];
            dmot6 = [zeros(1,6); diff(mot6)];
            mot24 = [mot6, dmot6, mot6.^2, dmot6.^2];
            dct_cols = confounds_table.Properties.VariableNames(startsWith(confounds_table.Properties.VariableNames, 'cosine'));
            
            nuisance = [mot24, table2array(confounds_table(:, [{'white_matter', 'csf'}, dct_cols]))];
            nuisance(isnan(nuisance)) = 0;
            X = [ones(size(nuisance, 1), 1), nuisance];
            
            % Regression and variance check
            residuals = bold_data - X * (X \ bold_data);
            v_std = std(residuals);
            good_vox = v_std > 0 & ~isnan(v_std);
            A = nan(size(residuals'));
            A(good_vox, :) = normalize(residuals(:, good_vox), 1)';
            
            ts_sub = A(subcort_all_idx, :); 
            ts_ctx = A(ctx_vertex_idx, :);
            
            % Temporal Correlation
            FC_vox_profile_all = nan(size(ts_sub, 1), size(ts_ctx, 1));
            FC_vox_profile_all(good_vox(subcort_all_idx), good_vox(ctx_vertex_idx)) = ...
                corr(ts_sub(good_vox(subcort_all_idx), :)', ts_ctx(good_vox(ctx_vertex_idx), :)', 'Rows', 'complete');
            
            FC_vox_profile = atanh(FC_vox_profile_all(local_idx_into_subcort, :));
            n_vox = size(FC_vox_profile, 1);
            corr_vals_left = zeros(n_vox, 2); corr_vals_right = zeros(n_vox, 2);
            
            for i = 1:n_vox
                row_prof = FC_vox_profile(i, :);
                
                corr_vals_left(i,1)  = corr(row_prof', ref_hip_l', 'Rows', 'complete');
                corr_vals_left(i,2)  = corr(row_prof', ref_amyg_l', 'Rows', 'complete');
                corr_vals_right(i,1) = corr(row_prof', ref_hip_r', 'Rows', 'complete');
                corr_vals_right(i,2) = corr(row_prof', ref_amyg_r', 'Rows', 'complete');
            end
    
            out_struct = struct('subj_id', subj_id, 'ses_name', ses_name, 'run', run, ...
                'voxel_indices', amyg_hip_vox_idx, 'subcortical_labels', parc(amyg_hip_vox_idx), ...
                'corr_left', double(corr_vals_left), 'corr_right', double(corr_vals_right));
            
            save(fullfile(output_dir, sprintf('CorrVals_AmygHipOnly_sub-%s_%s_run-%s.mat', subj_id, ses_name, run)), '-struct', 'out_struct');
        end
    end
end
output_table = cell2table(dropout_log(2:end,:), 'VariableNames', dropout_log(1,:));
writetable(output_table, fullfile(output_dir, 'dropout_log.csv'));