spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clear matlabbatch;

load('regressors_GLM_0405');
% % sub-10 제외.
% sbj_id_list(sbj_id_list==10)=[];
%sub-07 추가
sbj_id_list=sort([sbj_id_list 7]);


TYPE='FS60'; % 'HBT' or 'FS60'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\' TYPE];

beta_sbj1='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_glm\240426_anat\ODT\1st_Level\sub-01\beta_0001.nii';

% %% ODT
% for sbj_i=2:numel(sbj_id_list)
%     sbj_n = sbj_id_list(sbj_i);
%     beta_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_glm\240426_anat\ODT\1st_Level\' sprintf('sub-%.2d',sbj_n)];
%     betas=dir(fullfile(beta_path,'beta_*.nii'));
%     %% betas
%     for b=1:height(betas)
%         curr_beta = fullfile(betas(1).folder,betas(b).name);
%         clear matlabbatch;
%         matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {strcat(beta_sbj1,',1')};
%         matlabbatch{1}.spm.spatial.coreg.estwrite.source = {strcat(curr_beta,',1')};
%         matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
%         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
%         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
%         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%         matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
%         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
%         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
%         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
%         matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'coreg_';
% 
%         spm_jobman('run', matlabbatch)
% 
%     end
% end
% 
% 

%% mask
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'*.nii'));
    for m=1:height(local_masks)
        curr_mask=fullfile(local_masks(1).folder,local_masks(m).name);
        clear matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {strcat(beta_sbj1,',1')}; % 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\data_fmri_glm\240403\main\obj_show\1st_Level\sub-01\beta_0001.nii,1'};
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