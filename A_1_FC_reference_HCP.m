clear; clc;

Data = '/data/datasets/connectome/WU-Minn_HCP/3T_rfMRI_REST_fix';
MotionRoot = '/data/datasets/connectome/WU-Minn_HCP/RS-fMRI1_preproc';
Top = '/data/project';
Path = fullfile(Top, 'HCPDfMRI/fmriresults01');
Sub = dir(fullfile(Data, '??????'));

addpath(genpath(Path));
addpath(genpath(fullfile(Top, 'tools')));

cifti_label = cifti_read(fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii'));
parc = cifti_label.cdata;
template_scalar = cifti_label;
template_scalar.diminfo{2} = cifti_diminfo_make_scalars(1);

subcort_idx = find(parc >= 1 & parc <= 16);
cortex_idx = find(parc >= 17);

output_dir = fullfile(Top, 'derivatives/FC_reference_HCP');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

TR = 0.72;
hp_cut = 0.008; 
trim = 4;

%% Scan 1

FC_all = zeros(size(subcort_idx,1),size(cortex_idx,1));
N_used = 0;

for s = 1:height(Sub)
    subj_idx = s;
    sid = Sub(subj_idx).name;
    
    runs = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL'};
    clean_ts = {};
    
    for r = 1:2
        run = runs{r};

        ts_path = fullfile(Data, sid, 'MNINonLinear/Results', run, [run '_Atlas_MSMAll_hp0_clean.dtseries.nii']);
        
        ts_mov = fullfile(MotionRoot, sid, 'MNINonLinear/Results', run);
        fd_path = fullfile(ts_mov, 'Movement_RelativeRMS.txt');
        mot_path = fullfile(ts_mov, 'Movement_Regressors.txt');
       
        fd = load(fd_path);
        if (sum(fd > 0.5) / numel(fd)) > 0.2, continue; end
        
        img = cifti_read(ts_path); 
        T = img.cdata';
        bad = find(fd > 0.5);
        frames = unique([bad-1; bad; bad+1; bad+2]); 
        frames(frames < 1 | frames > size(T,1)) = [];
        if ~isempty(frames)
            good = setdiff(1:size(T,1), frames);
            T(frames,:) = interp1(good, T(good,:), frames, 'linear', 'extrap');
        end
        
        T = T(trim+1:end, :);
        mot = load(mot_path); mot = mot(trim+1:end, :);
        Nt = size(T,1);
        K = floor(2 * Nt * TR * hp_cut);
        DHP = cos(pi * ((2*(0:Nt-1)' + 1) .* (1:K)) / (2*Nt));
        X = [ones(Nt,1), mot, mot.^2, DHP];
        
        resid = T - X * (X \ T);
        clean_ts{end+1} = normalize(resid, 1);
    end
    

    if numel(clean_ts) == 2
        A = [clean_ts{1}; clean_ts{2}];
        
        ts_sub = A(:, subcort_idx);
        ts_ctx = A(:, cortex_idx);
        
        r = corr(ts_sub, ts_ctx); 
        z = atanh(max(min(r, 1-1e-7), -1+1e-7));
        
        FC_all = FC_all + z;
        N_used = N_used + 1;
    else
        fprintf('S1 Subject %d excluded (QC/Missing): %s\n', subj_idx, sid);
    end
end

if N_used > 0
    FC_mean = FC_all ./ N_used;
    save(fullfile(output_dir, 'FC_reference_HCP_scan1.mat'), 'FC_mean', 'N_used');
end


%% Scan 2
FC_all = zeros(size(subcort_idx,1),size(cortex_idx,1));
N_used = 0;

for s = 1:height(Sub)
    subj_idx = s;
    sid = Sub(subj_idx).name;
    
    runs = {'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
    clean_ts = {};
    
    for r = 1:2
        run = runs{r};
        ts_path = fullfile(Data, sid, 'MNINonLinear/Results', run, [run '_Atlas_MSMAll_hp0_clean.dtseries.nii']);
        
        ts_mov = fullfile(MotionRoot, sid, 'MNINonLinear/Results', run);
        fd_path = fullfile(ts_mov, 'Movement_RelativeRMS.txt');
        mot_path = fullfile(ts_mov, 'Movement_Regressors.txt');
        
        
        fd = load(fd_path);
        if (sum(fd > 0.5) / numel(fd)) > 0.2, continue; end 
        
        img = cifti_read(ts_path); 
        T = img.cdata'; 
        
        bad = find(fd > 0.5);
        frames = unique([bad-1; bad; bad+1; bad+2]); 
        frames(frames < 1 | frames > size(T,1)) = [];
        
        if ~isempty(frames)
            good = setdiff(1:size(T,1), frames);
            T(frames,:) = interp1(good, T(good,:), frames, 'linear', 'extrap');
        end
        
        T = T(trim+1:end, :);
        mot = load(mot_path); mot = mot(trim+1:end, :);
        
        Nt = size(T,1);
        K = floor(2 * Nt * TR * hp_cut);
        DHP = cos(pi * ((2*(0:Nt-1)' + 1) .* (1:K)) / (2*Nt));
        X = [ones(Nt,1), mot, mot.^2, DHP];
        
        resid = T - X * (X \ T);
        clean_ts{end+1} = normalize(resid, 1);
    end
    
    if numel(clean_ts) == 2
        A = [clean_ts{1}; clean_ts{2}];
        
        ts_sub = A(:, subcort_idx);
        ts_ctx = A(:, cortex_idx);
        
        r = corr(ts_sub, ts_ctx); 
        z = atanh(max(min(r, 1-1e-7), -1+1e-7));
        
        FC_all = FC_all + z;
        N_used = N_used + 1;
    else
        fprintf('S2 Subject %d excluded (QC/Missing): %s\n', subj_idx, sid);
    end
end

if N_used > 0
    FC_mean = FC_all ./ N_used;
    save(fullfile(output_dir, 'FC_reference_HCP_scan2.mat'), 'FC_mean', 'N_used');
end