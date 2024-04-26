%-----------------------------------------------------------------------
% Job saved on 25-Apr-2024 15:12:03 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clear matlabbatch;



load('regressors_GLM_0405');

TYPE='HBT'; % 'HBT' or 'FS60'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\' TYPE];


for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'*.nii'));

    beta_1=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\data_fmri_glm\240403\main\obj_show\1st_Level\' sprintf('sub-%.2d',sbj_n) '\beta_0001.nii'];

        for m=1:height(local_masks)
            curr_mask=fullfile(local_masks(1).folder,local_masks(m).name);
clear matlabbatch;
                matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {strcat(beta_1,',1')}; % 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\data_fmri_glm\240403\main\obj_show\1st_Level\sub-01\beta_0001.nii,1'};
                matlabbatch{1}.spm.spatial.coreg.estwrite.source = {strcat(curr_mask,',1')}; %{'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\HBT\s-01\lh_BODY.nii,1'};
                matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'coreg_';

                spm_jobman('run', matlabbatch)
        end
  
end