clear; clc;

Top = '/data/project';
Path = fullfile(Top, '/PNC/DataProcessed/sourcedata/freesurfer');
SubDirs = dir(fullfile(Path, 'sub_*'));
n_subj = length(SubDirs);

subject_ids = strings(n_subj,1);
icv_values = nan(n_subj,1);

for i = 1:n_subj
    subj_dir = SubDirs(i).name;
    subj_id = subj_dir;
    subject_ids(i) = subj_id;

    aseg_path = fullfile(Path, subj_dir, 'stats', 'aseg.stats');
    if ~isfile(aseg_path)
        warning('Missing aseg.stats for %s', subj_dir);
        continue;
    end

    fid = fopen(aseg_path, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    match_line = find(contains(lines{1}, 'EstimatedTotalIntraCranialVol'), 1);
    if ~isempty(match_line)
        parts = strsplit(lines{1}{match_line}, ',');
        icv_val = str2double(strtrim(parts{4}));
        icv_values(i) = icv_val;
    end
end

icv_table = table(subject_ids, icv_values, ...
    'VariableNames', {'src_subject_id','ICV'});

writetable(icv_table, fullfile(Top, 'derivatives/qc_metrics/icv_values_PNC.csv'));
