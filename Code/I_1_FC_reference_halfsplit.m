
clear; clc;


Top_Path = '/data/project';
Tools_Path = fullfile(Top_Path, 'tools');
addpath(genpath(Tools_Path));

temp_output_dir = fullfile(Top_Path, 'derivatives', 'temp_reliability_files');
output_dir = fullfile(Top_Path, 'derivatives', 'FC_reference_HCP');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

N_iterations = 100;

all_files = dir(fullfile(temp_output_dir, '*_Zmap_optimized.mat'));


first_data = load(fullfile(temp_output_dir, all_files(1).name));
Z_dims = size(first_data.Z);

Z_total = zeros(Z_dims);
for k = 1:numel(all_files)
    data = load(fullfile(temp_output_dir, all_files(k).name));
    Z_total = Z_total + data.Z;
end

N_subjects = numel(all_files);

Atlas_Path = fullfile(Top_Path, 'Atlas');
atlas_file = fullfile(Atlas_Path, 'Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
cifti_label = cifti_read(atlas_file); parc = cifti_label.cdata;
label_R_hip=1; label_R_amyg=2; label_L_hip=9; label_L_amyg=10;
amyghip_labels=[label_R_hip,label_R_amyg,label_L_hip,label_L_amyg];
amyghip_idx = find(ismember(parc, amyghip_labels));
parc_amyghip = parc(amyghip_idx);
mask_L_hip  = ismember(parc_amyghip, label_L_hip); mask_L_amyg = ismember(parc_amyghip, label_L_amyg);
mask_R_hip  = ismember(parc_amyghip, label_R_hip); mask_R_amyg = ismember(parc_amyghip, label_R_amyg);

reliability_coeffs = zeros(N_iterations, 4);

for i = 1:N_iterations
    shuffled_idx = randperm(N_subjects);
    split_point = floor(N_subjects / 2);
    groupA_idx = shuffled_idx(1:split_point);
    groupB_idx = shuffled_idx(split_point+1:end);
    
    sumZ_A = zeros(Z_dims);
    for k = 1:numel(groupA_idx)
        data = load(fullfile(temp_output_dir, all_files(groupA_idx(k)).name));
        sumZ_A = sumZ_A + data.Z;
    end
    
    sumZ_B = Z_total - sumZ_A;
    
    meanZ_A = sumZ_A / numel(groupA_idx);
    meanZ_B = sumZ_B / numel(groupB_idx);
    
    get_profile = @(Zmat, mask) tanh(mean(Zmat(mask, :), 1, 'omitnan'));
    profile_A_Lhip = get_profile(meanZ_A, mask_L_hip); profile_A_Lamyg = get_profile(meanZ_A, mask_L_amyg);
    profile_A_Rhip = get_profile(meanZ_A, mask_R_hip); profile_A_Ramyg = get_profile(meanZ_A, mask_R_amyg);
    profile_B_Lhip = get_profile(meanZ_B, mask_L_hip); profile_B_Lamyg = get_profile(meanZ_B, mask_L_amyg);
    profile_B_Rhip = get_profile(meanZ_B, mask_R_hip); profile_B_Ramyg = get_profile(meanZ_B, mask_R_amyg);
    
    reliability_coeffs(i, 1) = corr(profile_A_Lhip(:), profile_B_Lhip(:));
    reliability_coeffs(i, 2) = corr(profile_A_Lamyg(:), profile_B_Lamyg(:));
    reliability_coeffs(i, 3) = corr(profile_A_Rhip(:), profile_B_Rhip(:));
    reliability_coeffs(i, 4) = corr(profile_A_Ramyg(:), profile_B_Ramyg(:));
    
end

mean_reliability = mean(reliability_coeffs, 1); std_reliability = std(reliability_coeffs, 0, 1);

save(fullfile(output_dir, 'FC_reference_HCP_scan1_SplitHalfReliability.mat'), ...
     'reliability_coeffs', 'mean_reliability', 'std_reliability', 'N_iterations');