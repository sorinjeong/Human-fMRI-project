% -- Modified by Sorin Jeong, June 2024

% Result Report
% Random Effect Analysis (2nd-level Model for Between-Subject Analysis
clear; clc; close all;
%set path
addpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');
root_path = 'C:\Users\User\Desktop\JSR';%90번컴
addpath(genpath(root_path));
rehash path

%% Input!!

% %% phase of Interest
% % main: 1) 'obj_show', 2) 'choice', 3) 'moving'
% % ODT: 'ODT'
% 
% phase = 'obj_show';

%% 1) Defining pathway and subjects
% load regressors
load('reg_new_glm_0619.mat');
regressor=new_reg;
% set-up directory
path_in_1st = fullfile(root_path,'new_ocat_bids', 'derivatives'); % fmriprep's output file (before smoothing)
path_out_1st = fullfile(root_path, 'new_ocat_glm/main/1st_Level');

if ~exist(path_out_1st,"dir"); mkdir(path_out_1st);end
addpath(path_out_1st)

% create subjects beta and contrast result
beta_contrast_table = cell(length(sbj_id_list), 3);

%% 2) Initialize SPM

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;

%% 3) 
    for i = 1:length(sbj_id_list)
        c_sbj = sprintf('sub-%.2d', sbj_id_list(i)); disp(c_sbj)

        sbj_dir = fullfile(path_in_1st, c_sbj);
        current_beta_out = fullfile(path_out_1st,c_sbj);mkdir(current_beta_out);

        %% 4) 1st GLM - Extract Scan
        % read 4D-scan
        file_in = dir(fullfile(sbj_dir,'func',[c_sbj '*preproc_bold.nii.gz']));
        nii_info = niftiinfo(fullfile(file_in.folder,file_in.name));
        nii_data = niftiread(nii_info);

        % trim scan, extract only main task phase
            scan_start = regressor{1,i}.scan_num(1);
            scan_end = regressor{1,i}.scan_num(2);
  
            scan_range = scan_start:scan_end;
 
        if scan_range(end)>nii_info.ImageSize(end)
            scan_range(end) = nii_info.ImageSize(end);
        end
        extracted_data = nii_data(:, :, :, scan_range);


        % create nifti file
        nii_info.ImageSize = size(extracted_data);
        current_nifti = fullfile(path_out_1st,[c_sbj '_main_trimmed.nii']);
        niftiwrite(extracted_data, current_nifti, nii_info);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 5) 1st GLM - Make multiple regressor
        % event regressors
        valid_idx=~cellfun(@isempty, regressor{1,i}.onset);
        names = regressor{1,i}.name(valid_idx);
        onsets = regressor{1,i}.onset(valid_idx);
        durations = num2cell(zeros(size(regressor{1,i}.duration(valid_idx))));
        save(fullfile(path_out_1st, ['../', c_sbj, '_main_regressor_240619.mat']), 'names', 'onsets', 'durations');



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % movement regressors
        movereg_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        conf_file = dir(fullfile(path_in_1st, c_sbj,'func', ['*task-OCAT', '*confounds*.tsv']));   % find the file
        data1 = readtable(fullfile(conf_file.folder, conf_file.name), 'FileType', 'text', 'Delimiter', '\t');

        trans_x = find(strcmp(data1.Properties.VariableNames, movereg_names{1}));
        trans_y = find(strcmp(data1.Properties.VariableNames, movereg_names{2}));
        trans_z = find(strcmp(data1.Properties.VariableNames, movereg_names{3}));
        rot_x = find(strcmp(data1.Properties.VariableNames, movereg_names{4}));
        rot_y = find(strcmp(data1.Properties.VariableNames, movereg_names{5}));
        rot_z = find(strcmp(data1.Properties.VariableNames, movereg_names{6}));

        % remove the first row, 26-31 columns --> movement regressors    R_mov1 = fillmissing(R_mov1, ?nearest?);    R_mov1(~isfinite(R_mov1))=0;
        current_mov_reg = data1{scan_range, [trans_x,trans_y,trans_z,rot_x,rot_y,rot_z]};
        current_mov_reg = fillmissing(current_mov_reg, 'nearest');
        current_mov_reg(~isfinite(current_mov_reg)) = 0;

        R = current_mov_reg;
        R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};

        save(fullfile(path_out_1st, ['../', c_sbj, '_movements.mat']), 'R', 'R_names');

        %% 6) 1st GLM - fMRI model specification
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_spec.dir = {current_beta_out};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {current_nifti};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(path_out_1st, ['../', c_sbj, '_main_regressor_240619.mat'])};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(path_out_1st, ['../', c_sbj, '_movements.mat'])};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;


        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1; %Model Interactions (Volterra) / 1==no interaction , 2==interaction
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; %autoregrssive ; AR(1) is the default option (by AHN)
        spm_jobman('run', matlabbatch)

        %% 7) 1st GLM - Model estimation
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(current_beta_out, 'SPM.mat') };
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', matlabbatch)

        %% 8) 1st GLM - Contrast manager
        % matlabbatch = [];
        % matlabbatch{1}.spm.stats.con.spmmat  = {fullfile(current_beta_out, 'SPM.mat') };
        % for c = 1:length(contrast{1,2})
        %     matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = contrast{1,2}{c};
        %     matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = contrast{2,2}{c,:};
        %     matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
        % 
        % end
        % matlabbatch{1}.spm.stats.con.delete = 1;
        % 
        % spm_jobman('run', matlabbatch)


    end




