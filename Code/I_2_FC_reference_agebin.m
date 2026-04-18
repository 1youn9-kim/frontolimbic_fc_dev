clear; clc;
Top = '/data/projects';
Data = '/data/datasets/connectome/WU-Minn_HCP/3T_rfMRI_REST_fix';
MotionRoot = '/data/datasets/connectome/WU-Minn_HCP/RS-fMRI1_preproc';
addpath(genpath(fullfile(Top, 'tools')));
output_dir = fullfile(Top, 'derivatives/FC_reference_HCP');
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end
pheno = readtable(fullfile(Top, 'HCPYA_pheno.csv'));
Sub = dir(Data);
Sub = Sub([Sub.isdir] & ~ismember({Sub.name}, {'.', '..'}));
atlas_path = fullfile(Top, 'Atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
cifti_label = cifti_read(atlas_path);
parc = cifti_label.cdata;
medial = [42,57,59,73,74,75,76,77,78,79,80,81,85,86,88,104,106,109,180,181,182,195,196,222,237,239,253,254,255,256,257,258,259,260,261,265,266,268,284,286,289,360,361,362,375,376];
lateral = [26,27,28,60,84,87,93,94,96,97,98,100,101,102,103,105,107,112,113,186,187,206,207,208,240,264,267,273,274,277,278,280,281,282,283,285,287,292,293,366,367,82,83,89,90,91,92,95,99,108,110,114,262,263,269,270,271,272,275,276,279,288,290,294];
ctx_vertex_idx = find(ismember(parc, [medial, lateral]));
subcort_all_idx = find(parc >= 1 & parc <= 16);
TR = 0.72; hp_cut = 0.008; trim = 4;
bin_names = {'22_26', '26_30', '30_35'};
for b = 1:3
    if b == 1
        bin_idx = find(pheno.Age_in_Yrs >= 22 & pheno.Age_in_Yrs < 26);
    elseif b == 2
        bin_idx = find(pheno.Age_in_Yrs >= 26 & pheno.Age_in_Yrs < 30);
    else
        bin_idx = find(pheno.Age_in_Yrs >= 30 & pheno.Age_in_Yrs <= 35);
    end
    
    bin_subs = pheno.Subject(bin_idx);
    FC_sum = zeros(length(subcort_all_idx), length(ctx_vertex_idx));
    valid_count = 0;
    
    for s = 1:length(bin_subs)
        sid = num2str(bin_subs(s));
        if ~exist(fullfile(Data, sid), 'dir')
            continue;
        end
        
        runs = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL'};
        clean_ts = {};
        
        for r = 1:2
            run = runs{r};
            ts_path = fullfile(Data, sid, 'MNINonLinear/Results', run, [run '_Atlas_MSMAll_hp2000_clean.dtseries.nii']);
            ts_mov = fullfile(MotionRoot, sid, 'MNINonLinear/Results', run);
            
            if ~exist(ts_path, 'file')
                continue; 
            end
            
            fd = load(fullfile(ts_mov, 'Movement_RelativeRMS.txt'));
            if (sum(fd > 0.5) / numel(fd)) > 0.2
                continue; 
            end
            
            T = cifti_read(ts_path).cdata';
            bad = find(fd > 0.5);
            frames = unique([bad-1; bad; bad+1; bad+2]); 
            frames(frames < 1 | frames > size(T,1)) = [];
            
            if ~isempty(frames)
                good = setdiff(1:size(T,1), frames);
                T(frames,:) = interp1(good, T(good,:), frames, 'linear', 'extrap');
            end
            
            T = T(trim+1:end, :);
            mot = load(fullfile(ts_mov, 'Movement_Regressors.txt')); 
            mot = mot(trim+1:end, :);
            
            Nt = size(T,1);
            K = floor(2 * Nt * TR * hp_cut);
            DHP = cos(pi * ((2*(0:Nt-1)' + 1) .* (1:K)) / (2*Nt));
            X = [ones(Nt,1), mot, mot.^2, DHP];
            
            resid = T - X * (X \ T);
            clean_ts{end+1} = normalize(resid, 1);
        end
        
        if numel(clean_ts) == 2
            A = [clean_ts{1}; clean_ts{2}]';
            ts_sub = A(subcort_all_idx, :);
            ts_ctx = A(ctx_vertex_idx, :);
            FC_vox = atanh(corr(ts_sub', ts_ctx'));
            FC_sum = FC_sum + FC_vox;
            valid_count = valid_count + 1;
        end
    end
    
    if valid_count > 0
        FC_mean = FC_sum / valid_count;
        out_name = fullfile(output_dir, sprintf('FC_reference_HCP_scan1_%s.mat', bin_names{b}));
        save(out_name, 'FC_mean');
    end
end
%% Visualiazation
amyg_sub_idx = ismember(parc(subcort_all_idx), [2, 10]);
hipp_sub_idx  = ismember(parc(subcort_all_idx), [1, 9]);

amyg_all = zeros(length(bin_names), length(ctx_vertex_idx));
hipp_all = zeros(length(bin_names), length(ctx_vertex_idx));

for b = 1:length(bin_names)
    load(fullfile(output_dir, sprintf('FC_reference_HCP_scan1_%s.mat', bin_names{b})), 'FC_mean');
    
    amyg_fc = mean(FC_mean(amyg_sub_idx, :), 1);
    hipp_fc  = mean(FC_mean(hipp_sub_idx, :), 1);
    
    amyg_all(b, :) = amyg_fc;
    hipp_all(b, :) = hipp_fc;
    
    out_map_amyg = zeros(size(parc, 1), 1);
    out_map_amyg(ctx_vertex_idx) = zscore(amyg_fc, 0, 2);
    
    out_cifti_amyg = cifti_label;
    out_cifti_amyg.cdata = out_map_amyg;
    out_cifti_amyg.diminfo{2}.type = 'scalars';
    if isfield(out_cifti_amyg.diminfo{2}, 'labels')
        out_cifti_amyg.diminfo{2} = rmfield(out_cifti_amyg.diminfo{2}, 'labels');
    end
    cifti_write(out_cifti_amyg, fullfile(output_dir, sprintf('FC_reference_Amyg_%s.dscalar.nii', bin_names{b})));
    
    out_map_hipp = zeros(size(parc, 1), 1);
    out_map_hipp(ctx_vertex_idx) = zscore(hipp_fc, 0, 2);
    
    out_cifti_hipp = cifti_label;
    out_cifti_hipp.cdata = out_map_hipp;
    out_cifti_hipp.diminfo{2}.type = 'scalars';
    if isfield(out_cifti_hipp.diminfo{2}, 'labels')
        out_cifti_hipp.diminfo{2} = rmfield(out_cifti_hipp.diminfo{2}, 'labels');
    end
    cifti_write(out_cifti_hipp, fullfile(output_dir, sprintf('FC_reference_Hipp_%s.dscalar.nii', bin_names{b})));
end

amyg_corr = corr(amyg_all');
hipp_corr = corr(hipp_all');

writematrix(amyg_corr, fullfile(output_dir, 'FC_reference_Amyg_crosscorr.csv'));
writematrix(hipp_corr, fullfile(output_dir, 'FC_reference_Hipp_crosscorr.csv'));