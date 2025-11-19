clear; clc;

Top = '/data/project';
Path = fullfile(Top, 'HCPDfMRI/fmriresults01');
addpath(genpath(fullfile(Top, 'tools')));

subs = readtable(fullfile(Top, 'HCPDfMRI/ndar_subject01.txt'));
n_subs = height(subs);

results = cell(n_subs, 2);

for s = 1:n_subs
    subj_id = subs.src_subject_id{s};
    
    sub_dir_name = sprintf('%s_V1_MR', subj_id);
    aseg_path = fullfile(Path, sub_dir_name, 'T1w', subj_id, 'stats', 'aseg.stats');
    
    icv_val = NaN;
    
    fid = fopen(aseg_path, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    
    match_line = find(contains(lines{1}, 'EstimatedTotalIntraCranialVol'), 1);
    
    if ~isempty(match_line)
        parts = strsplit(lines{1}{match_line}, ',');
        icv_val = str2double(strtrim(parts{4}));
    end

    
    results{s, 1} = subj_id;
    results{s, 2} = icv_val;
end

icv_table = cell2table(results, 'VariableNames', {'src_subject_id', 'ICV'});
writetable(icv_table, fullfile(Top, 'derivatives/qc_metrics', 'icv_values_HCPD.csv'));