clear; clc;
Top = '/data/projects/punim2400'; 
bcp_project_dir = fullfile(Top, 'BABY/image03/BCP');
Path = fullfile(bcp_project_dir, 'derivatives', 'nibabies', 'sourcedata', 'freesurfer');
subject_list_file = fullfile(bcp_project_dir, 'final_sublist.txt');

subjects_table = readtable(subject_list_file, 'ReadVariableNames', false);
n_subj = height(subjects_table);
results = {};

for i = 1:n_subj
    subj_id = sprintf('%06d', subjects_table.Var1(i));
    subj_fs_dirs = dir(fullfile(Path, ['sub-' subj_id '*']));
    subj_fs_dirs = subj_fs_dirs([subj_fs_dirs.isdir]);
    
    if isempty(subj_fs_dirs)
        results(end+1, :) = {['sub-' subj_id], 'N/A', NaN};
        continue;
    end
    
    for j = 1:numel(subj_fs_dirs)
        ses_dir_name = subj_fs_dirs(j).name;
        ses_tok = regexp(ses_dir_name, '(ses-[A-Za-z0-9]+)', 'tokens', 'once');
        if ~isempty(ses_tok), ses_name = ses_tok{1}; else, ses_name = 'N/A'; end
        
        brainvol_path = fullfile(Path, ses_dir_name, 'stats', 'brainvol.stats');
        aseg_path = fullfile(Path, ses_dir_name, 'stats', 'aseg.stats');
        
        icv_val = NaN;
        
        if isfile(brainvol_path)
            fid = fopen(brainvol_path, 'r');
            if fid ~= -1
                lines = textscan(fid, '%s', 'Delimiter', '\n');
                fclose(fid);
                idx = find(contains(lines{1}, 'MaskVol,'), 1);
                if ~isempty(idx)
                    parts = strsplit(lines{1}{idx}, ',');
                    if length(parts) >= 4
                        icv_val = str2double(strtrim(parts{4}));
                    end
                end
            end
        end

        if isnan(icv_val) && isfile(aseg_path)
            fid = fopen(aseg_path, 'r');
            if fid ~= -1
                lines = textscan(fid, '%s', 'Delimiter', '\n');
                fclose(fid);
                idx = find(contains(lines{1}, 'EstimatedTotalIntraCranialVol'), 1);
                if ~isempty(idx)
                    parts = strsplit(lines{1}{idx}, ',');
                    if length(parts) >= 4
                        icv_val = str2double(strtrim(parts{4}));
                    end
                end
            end
        end
        
        results(end+1, :) = {['sub-' subj_id], ses_name, icv_val};
    end
end

icv_table = cell2table(results, 'VariableNames', {'src_subject_id', 'session', 'ICV'});
icv_table.src_subject_id = string(icv_table.src_subject_id);
icv_table.session = string(icv_table.session);

output_dir = fullfile(Top, 'derivatives', 'BCP_FC_Analysis');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
writetable(icv_table, fullfile(output_dir, 'icv_values_BCP.csv'));