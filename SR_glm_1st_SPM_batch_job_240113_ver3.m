clear; clc; close all;


%% Initialise SPM defaults
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%-----------------------------------------------------------------------
% Job saved on 13-Jan-2024 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
path_in = 'D:\fMRI\OCAT_DIR\data\data_fmri_bids\derivatives'; % fmriprep's output file (before smoothing)
path_out = 'D:\fMRI\OCAT_DIR\data\lev-1st_ver4'; % save betas this folder
addpath(path_in);
condition_path = 'D:\fMRI\OCAT_DIR\data\data_for_glm';addpath(condition_path);
load('regressor_all_all'); 


%% 
for i = 1:length(sbj_id_list)
    c_sbj = sprintf('sub-%.2d', sbj_id_list(i));
    
    sbj_dir = fullfile(path_in, c_sbj);
    current_beta_out = fullfile(path_out,c_sbj);mkdir(current_beta_out);

%% Version : smoothing O --> if you want to add a smoothing step, run the code starting from this block

%     file_in = dir(fullfile(sbj_dir,'func',[c_sbj '*preproc_bold.nii']));
%     current_file_in = cellstr(fullfile(file_in.folder,file_in.name));
% 
% 
% matlabbatch{1}.spm.spatial.smooth.data = {current_file_in};
% matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
% matlabbatch{1}.spm.spatial.smooth.dtype = 0;
% matlabbatch{1}.spm.spatial.smooth.im = 0;
% matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% 
% 

%% Version : smoothing X

    file_in = dir(fullfile(sbj_dir,'func',['s' c_sbj '*preproc_bold.nii']));
    current_file_in = cellstr(fullfile(file_in.folder,file_in.name));
 
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_spec.dir = {current_beta_out};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = current_file_in;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(condition_path,[c_sbj, '_multiple_conditions.mat'])};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(condition_path,[c_sbj, '_movereg.mat'])};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1; %Model Interactions (Volterra) / 1==no interaction , 2==interaction 
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
spm_jobman('run', matlabbatch)

%%
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_est.spmmat = { fullfile(current_beta_out, 'SPM.mat') };
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch)

%%
matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat  = { fullfile(current_beta_out, 'SPM.mat') };
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'corr-incorr';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0 1 -1 1 -1 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'corr_YES-corr_NO';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'corr_NO-corr_YES';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 -1 0 1 0];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch)

end


%%%%% 2:37am : 지금 계속 #31번피험자 (28번째) 가 fMRI model specification 실패중임. 일단
%%%%% 제외하고 2nd 진행

%%%%% 2024/01/19 8pm : #31번 smoothing 새로 하고 (전에도 smoothing 다시 했었지만 그 땐
%%%%% 안됐었는데..) 한 블럭 한블럭씩 하니까 됨. smoothing 한번 돌리고 specification 돌리고 이런식으로.
