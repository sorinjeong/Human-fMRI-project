clear; clc
addpath('C:\Users\Leelab\Documents\MATLAB\spm12');
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clear matlabbatch;

load('regressors_GLM_0513');
ODT=odt_reg{1,1}.ODT;

TYPE='FS60'; % 'HBT' or 'FS60' or 'extra'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_bids\derivatives\HPC_seg\' TYPE];
beta_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_glm\240512_new\ODT\1st_Level';

%% ODT
for sbj_i=2:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    betas=dir(fullfile(beta_path,sprintf('sub-%.2d',sbj_n),'beta_*.nii'));

    mask_sbj=fullfile(mask_path,sprintf('s-%.2d',sbj_n));
    masks=dir(fullfile(mask_sbj,'*.nii'));
    %% betas
    roi_dir=fullfile(beta_path,'roi',sprintf('sub-%.2d',sbj_n));mkdir(roi_dir);

    for b=1:height(betas)-7 %movement regressors(6) and intercept(1)
        ref_beta=fullfile(betas(1).folder,betas(b).name);

        %% masks
        for m=1:height(masks)
            curr_mask = fullfile(masks(1).folder,masks(m).name);

            clear matlabbatch;
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {strcat(ref_beta,',1')};
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {strcat(curr_mask,',1')};
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = ['coreg_' ODT.regress_name{b} '_'];

            spm_jobman('run', matlabbatch)

        end
    end
end