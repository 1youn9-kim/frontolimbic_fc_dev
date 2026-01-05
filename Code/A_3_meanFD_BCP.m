clear; clc;

Top = '/data/projects/punim2400';
bcp_project_dir = fullfile(Top, 'BABY/image03/BCP');
search_root = fullfile(bcp_project_dir, 'derivatives', 'nibabies');
output_dir = fullfile(Top, 'derivatives', 'BCP_FC_Analysis');

dd = dir(fullfile(search_root, '**', '*_desc-confounds_timeseries.tsv'));
results = {};

for i = 1:length(dd)
    filepath = fullfile(dd(i).folder, dd(i).name);
    
    [~, filename, ~] = fileparts(filepath);
    
    subj_tok = regexp(filename, 'sub-(\d{1,6})', 'tokens', 'once');
    ses_tok  = regexp(filename, 'ses-([A-Za-z0-9]+)', 'tokens', 'once');
    run_tok  = regexp(filename, 'run-([0-9]+)', 'tokens', 'once');
    
    subj_id = ['sub-' subj_tok{1}];
    ses_name = ['ses-' ses_tok{1}];
    run_id   = string(run_tok{1});
    
    T = readtable(filepath, 'FileType', 'text', 'Delimiter', '\t', ...
            'TreatAsEmpty', {'n/a','NA'}, 'VariableNamingRule', 'preserve');

    fd_col = T.framewise_displacement;
    
    if isempty(fd_col) || height(T) < 2, continue; end
    if isnan(fd_col(1)), fd_col(1) = 0; end
    
    mean_fd_run = mean(fd_col, 'omitnan');
    
    results(end+1, :) = {subj_id, ses_name, run_id, mean_fd_run};
end

if ~isempty(results)
    G = cell2table(results, 'VariableNames', {'src_subject_id', 'session', 'run', 'MeanFD'});
    writetable(G, fullfile(output_dir, 'mean_fd_run_BCP.csv'));
end