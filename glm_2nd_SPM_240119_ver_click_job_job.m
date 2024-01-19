clear; clc; close all;


%% Initialise SPM defaults
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%-----------------------------------------------------------------------
% Job saved on 19-Jan-2024 20:30:17 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
path_in = 'D:\fMRI\OCAT_DIR\data\lev-1st_ver4'; % 1st level GLM - output file 
path_out = 'D:\fMRI\OCAT_DIR\data\lev-2nd-ver5'; % save betas this folder
addpath(path_in);
condition_path = 'D:\fMRI\OCAT_DIR\data\data_for_glm';addpath(condition_path);
load('regressor_all_all'); 

%% set conditions

cond_name = {'corr-incorr','corr_YES-corr_NO', 'corr_NO-corr_YES'};
cond_num = {'con_0001','con_0002','con_0003'};

for c=1:3
    fprintf('%s%s\n','current condition : ',cond_name{c});
    
%% set path
 current_path_out = fullfile(path_out,cond_name{c}); mkdir(current_path_out);
%% set files
current_scans = char(arrayfun(@(i) [path_in,'\sub-', sprintf('%02d', i), '\',cond_num{c},'.nii,1'], sbj_id_list, 'UniformOutput', false));

%% 
matlabbatch=[];
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(current_path_out);
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(string(current_scans));
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'corr-incorr';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch)

end