cd 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code'
addpath('C:\Users\User\Documents\MATLAB\');
load('regressors_GLM_0405');

TYPE='HBT'; % 'HBT' or 'FS60'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\' TYPE];


roi=struct;

for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'coreg_*.nii'));

    beta_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\data_fmri_glm\240403\main\obj_show\1st_Level\' sprintf('sub-%.2d',sbj_n)];
    betas=dir(fullfile(beta_path,'beta_*.nii'));

    for b=1:height(betas)
        curr_beta = fullfile(betas(1).folder,betas(b).name);
        for m=1:height(local_masks)

            curr_mask=niftiread(fullfile(local_masks(1).folder,local_masks(m).name));
            roi_coord = mask2coord(curr_mask);

            roi_beta = niftiread(curr_beta);
            roi_beta_2d = reshape(roi_beta,1001,300);

        roi.(sbj_n).hpc.roi_name_list(m)=local_masks(m).name;
        roi.(sbj_n).hpc.roi_pattern(m)=roi_beta_2d;

            
        end
    end
end



