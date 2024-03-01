clear; clc; close all;
%%% regressor 시간을 2초씩 뒤로 미뤘고, scan을 ODT때 양 끝을 제외한 scan 이용

%% Initialise SPM defaults
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%-----------------------------------------------------------------------
% Job saved on 13-Jan-2024 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
root_path = 'D:\fMRI\OCAT_DIR\data';

path_in = fullfile(root_path,'data_fmri_bids','derivatives'); % fmriprep's output file (before smoothing)
path_out = fullfile(root_path,'lev-1st_0125_ver2'); % save betas this folder
addpath(path_in);
condition_path = fullfile(root_path,'data_for_glm_240125');addpath(condition_path);
load('event_regressor_all2'); 


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
    main_task_info = niftiinfo(fullfile(file_in.folder,file_in.name));
    main_task_data = niftiread(main_task_info);

    %% extracting scan 
    main_start = behav_regressor_glm{1,4}{1,i}.scan_num(1);
    main_end = behav_regressor_glm{1,4}{1,i}.scan_num(2);
    extracted_data = main_task_data(:, :, :, main_start:main_end);
    
    % calculate main scan onset
    
    main_onset = main_task_info.PixelDimensions(4) * (main_start-1); % PixelDimensions(4) == TR을 의미
    
    current_nifti = main_task_info;
    current_nifti.ImageSize = size(extracted_data);

    niftiwrite(extracted_data, fullfile(path_out,'current_file_in.nii'), current_nifti);

    %% multiple regressor
    names = behav_regressor_glm{1, 4}{1, i}.regress_name;
    onsets = behav_regressor_glm{1, 4}{1, i}.regress_onset;
    durations = behav_regressor_glm{1, 4}{1, i}.regress_duration;

     save(fullfile(condition_path,[c_sbj, '_multiple_conditions.mat']), 'names', 'onsets', 'durations');
 
    % multiple condition
    movereg_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
    
    conf_file = dir(fullfile(path_in, c_sbj,'func', ['*task-OCAT', '*confounds*.tsv']));   % find the file
    [data1, header1, ] = tsvread(fullfile(conf_file.folder, conf_file.name));
    
    trans_x=strmatch(movereg_names{1},header1,'exact');
    trans_y=strmatch(movereg_names{2},header1,'exact');
    trans_z=strmatch(movereg_names{3},header1,'exact');
    rot_x=strmatch(movereg_names{4},header1, 'exact');
    rot_y=strmatch(movereg_names{5},header1,'exact');
    rot_z=strmatch(movereg_names{6},header1,'exact');
    
    % remove the first row, 26-31 columns --> movement regressors    R_mov1 = fillmissing(R_mov1, ?nearest?);    R_mov1(~isfinite(R_mov1))=0;
    current_mov_reg = data1(main_start+1:main_end+1, [trans_x,trans_y,trans_z,rot_x,rot_y,rot_z]);
    current_mov_reg = fillmissing(current_mov_reg, 'nearest');
    current_mov_reg(~isfinite(current_mov_reg)) = 0;
    
    %% 
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_spec.dir = {current_beta_out};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(fullfile(path_out,'current_file_in.nii'));
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(condition_path,[c_sbj, '_multiple_conditions.mat'])};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {current_mov_reg};
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
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1 1 -1 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
if sbj_id_list(i)==7
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 1 0];
end
if sbj_id_list(i)==25
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1 1 0];
end
if sbj_id_list(i)==28
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1 -1 0];
end
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'corr_YES-corr_NO';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'corr_NO-corr_YES';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [-1 0 1 0];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch)

end


%%%%% 2024/01/25 #7은 incorrect 없으므로 그냥 OCP면 모두 weight 1로 부여 

