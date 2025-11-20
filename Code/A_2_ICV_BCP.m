clear; clc;

Top = '/data/projects'; 
bcp_project_dir = fullfile(Top, 'BABY/image03/BCP');

Path = fullfile(bcp_project_dir, 'derivatives/sourcedata/freesurfer');
subject_list_file = fullfile(bcp_project_dir, 'final_sublist.txt');

subjects_table = readtable(subject_list_file, 'ReadVariableNames', false);
n_subj = height(subjects_table);
subjects_to_process = cell(n_subj, 1);
for i = 1:n_subj
    subjects_to_process{i} = sprintf('%06d', subjects_table.Var1(i));
end

results = {};

for i = 1:n_subj
    subj_id = subjects_to_process{i};
    
    subj_fs_dirs = dir(fullfile(Path, ['sub-' subj_id '_ses-*']));
    subj_fs_dirs = subj_fs_dirs([subj_fs_dirs.isdir]);
    
    for j = 1:numel(subj_fs_dirs)
        ses_dir_name = subj_fs_dirs(j).name;
        
        ses_tok = regexp(ses_dir_name, '(ses-[A-Za-z0-9]+)', 'tokens', 'once');
        if isempty(ses_tok), continue; end
        ses_name = ses_tok{1};
        
        aseg_path = fullfile(Path, ses_dir_name, 'stats', 'aseg.stats');
        
        if ~isfile(aseg_path), continue; end
        
        fid = fopen(aseg_path, 'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        match_line = find(contains(lines{1}, 'EstimatedTotalIntraCranialVol'), 1);
        
        icv_val = NaN;
        if ~isempty(match_line)
            parts = strsplit(lines{1}{match_line}, ',');
            icv_val = str2double(strtrim(parts{4}));
        end
        
        results(end+1, :) = {['sub-' subj_id], ses_name, icv_val};
    end
end

icv_table = cell2table(results, ...
    'VariableNames', {'src_subject_id', 'session', 'ICV'});

writetable(icv_table, fullfile(Top, 'derivatives', 'BCP_FC_Analysis', 'icv_values_BCP.csv'));