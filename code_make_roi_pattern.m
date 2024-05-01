clear;clc;

cd 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code'
addpath('C:\Users\User\Documents\MATLAB\');
load('regressors_GLM_0405');
% sub-10 제외.
sbj_id_list(sbj_id_list==10)=[];
%sub-07 추가
sbj_id_list=sort([sbj_id_list 7]);

TYPE='HBT'; % 'HBT' or 'FS60'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\' TYPE];


roi=struct;
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'coreg_*.nii'));

    beta_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_glm\240426_anat\ODT\1st_Level\' sprintf('sub-%.2d',sbj_n)];
    betas=dir(fullfile(beta_path,'coreg_beta_*.nii'));

    for b=1:height(betas)
        curr_beta = fullfile(betas(1).folder,betas(b).name);
        for m=1:height(local_masks)
            curr_mask=niftiread(fullfile(local_masks(1).folder,local_masks(m).name));
            roi_coord = mask2coord(curr_mask);

            roi_beta = niftiread(curr_beta);
            roi_beta_2d = reshape(roi_beta,416,354);

        roi.(sprintf('s%.2d',sbj_n)).hpc.(TYPE).roi_name_list{m}=local_masks(m).name;
        roi.(sprintf('s%.2d',sbj_n)).hpc.(TYPE).roi_pattern{m}=roi_beta_2d;

            
        end
    end
end


%% roi_beta reshape할 때 beta 마다 달라짐!!!

%ODT: 52*59*48 = 416*354 = 147264
%obj_show: 1001*300 = 300300



