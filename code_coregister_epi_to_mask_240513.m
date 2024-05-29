clear; clc
% addpath('C:\Users\Leelab\Documents\MATLAB\spm12');%69번컴은 leelab, 소린컴은 Leelab
addpath('C:\Users\User\Documents\MATLAB\spm12\');%90번
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clear matlabbatch;

load('regressors_GLM_0513');
ODT=odt_reg{1,1}.ODT;

TYPE='FS60'; % 'HBT' or 'FS60' or 'extra'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_bids\derivatives\HPC_seg\' TYPE];
beta_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_bids\derivatives';

%% ODT epi to mask coregister
        beta_names={'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2','post_forest_1','post_forest_2','post_city_1','post_city_2','pre_target','post_target'};


for sbj_i=22:numel(sbj_id_list)
% sbj_i=3;
    sbj_n = sbj_id_list(sbj_i);
    betas=dir(fullfile(beta_path,sprintf('sub-%.2d',sbj_n),'func','*_space-T1w_boldref.nii'));
    masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'*.nii'));

    %% betas
%     roi_dir=fullfile(beta_path,'roi',sprintf('sub-%.2d',sbj_n));mkdir(roi_dir);

    % for b=1:height(betas)-7 %movement regressors(6) and intercept(1)
        ref_beta=fullfile(betas(1).folder,betas(1).name);
        

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
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'coreg_';

            spm_jobman('run', matlabbatch)
            % 
            % %% move file
            % 
            % cor_mask=dir(fullfile(betas(1).folder,'coreg_*.nii'));
            % 
            % bn=beta_names(b);
            % roi_path=fullfile('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_glm\240512_new\ODT\1st_Level\roi_mask_new', TYPE, bn,sprintf('sub-%.2d',sbj_n));
            % mkdir(roi_path{:})
            % 
            % movefile(fullfile(cor_mask(end).folder,cor_mask(end).name),fullfile(roi_path{:},strcat(extractBetween(cor_mask(1).name,"coreg_","_beta"),"_",masks(m).name)))


        end
    
end

