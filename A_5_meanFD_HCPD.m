clear; clc;

Top = '/data/project';
Path = fullfile(Top, 'HCPDfMRI/fmriresults01');
Sub = dir(fullfile(Path, '*_MR'));
addpath(genpath(fullfile(Top, 'tools')));

subs = readtable(fullfile(Top, 'HCPDfMRI/ndar_subject01.txt'));
subs = sortrows(subs, "src_subject_id");

fd_vals = nan(height(subs), 1);
sub_ids = strings(height(subs), 1);

fd_out_csv = fullfile(Top, 'derivatives/qc_metrics/mean_fd_allsubs.csv');

for s = 1:height(subs)
    subj_id = subs.src_subject_id{s};
    subname = sprintf('%s_V1_MR', subj_id);
    sub_path = fullfile(Path, subname, 'MNINonLinear', 'Results');

    fd_files = dir(fullfile(sub_path, '*', 'Movement_RelativeRMS.txt'));

    fd_vals_each = nan(1, length(fd_files));
    for f = 1:length(fd_files)
        fd_series = readmatrix(fullfile(fd_files(f).folder, fd_files(f).name));
        fd_vals_each(f) = mean(fd_series);

    end

    sub_ids(s) = subj_id;
    fd_vals(s) = mean(fd_vals_each, 'omitnan');
end

T = table(sub_ids, fd_vals, 'VariableNames', {'src_subject_id', 'MeanFD'});
writetable(T, fd_out_csv);
