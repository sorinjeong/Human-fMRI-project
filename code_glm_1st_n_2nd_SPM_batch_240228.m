% -- Modified by Sorin Jeong, February 2024

% Result Report
% Random Effect Analysis (2nd-level Model for Between-Subject Analysis
clear; clc; close all;
%set path
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');
root_path = 'D:\fMRI\OCAT_DIR\data';
addpath(genpath(root_path));

%% Input!!
add_smoothing = 0;% if you want to add a smoothing step, put 1
run_1st_glm = 1; % if you want to run the 1st glm step, put 1
run_2nd_glm = 1; % if you want to run the 1st glm step, put 1
threshold_type = 'FWE' ;% 'FWE' or 'FDR' or 'none'
date = '240301';
% isROI = 1; % if you want to run the ROI analysis, put 1, or if not, put 0 for whole brain analysis

%% contrasts!!
contrast{1,1}='name: ';
contrast{1,2}={'correct';'incorrect';'Hit'; 'correct_rejection';'Corr-Incorr';'Incorr-Corr';'Hit-Corr_reject';'Corr_reject-Hit'};
contrast{2,1}='weight: ';
contrast{2,2}={'[0.5 0 0.5 0]'; '[0 0.5 0 0.5]'; '[1 0 0 0]'; '[0 0 1 0]'; '[0.5 -0.5 0.5 -0.5]'; '[-0.5 0.5 -0.5 0.5]'; '[1 0 -1 0]'; '[-1 0 1 0]'};
save(char(fullfile(root_path,date,'contrast_v1.mat')),"contrast");
%% 1) Defining pathway and subjects

% set-up directory
path_in_1st = fullfile(root_path,'data_fmri_bids','derivatives'); % fmriprep's output file (before smoothing)
path_out_1st = fullfile(root_path, date,'1st_Level'); % save betas this folder
path_out_2nd = fullfile(root_path,date,'2nd_Level'); % save 2nd lev glm results this folder
if ~exist(path_out_1st); mkdir(path_out_1st);end
if ~exist(path_out_2nd); mkdir(path_out_2nd);end

condition_path = fullfile(root_path,['data_for_glm_' date]);mkdir(condition_path);addpath(condition_path);

% load regressors
load('event_regressor_all2');

%% 2) Initialize SPM
spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;

%% 3) preprocessing - Spatial Smoothing
if run_1st_glm == 1
    for i = 1:length(sbj_id_list)
        c_sbj = sprintf('sub-%.2d', sbj_id_list(i));
        
        sbj_dir = fullfile(path_in_1st, c_sbj);
        current_beta_out = fullfile(path_out_1st,c_sbj);mkdir(current_beta_out);
        
        % yes_add_smoothing
        if add_smoothing == 1
            file_in = dir(fullfile(sbj_dir,'func',[c_sbj '*preproc_bold.nii']));
            current_file_in = cellstr(fullfile(file_in.folder,file_in.name));
            
            matlabbatch{1}.spm.spatial.smooth.data = {current_file_in};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            
            spm_jobman('run', matlabbatch)
        end
        
        %% 4) 1st GLM - Extract Scan
        % read 4D-scan
        file_in = dir(fullfile(sbj_dir,'func',['s' c_sbj '*preproc_bold.nii']));
        main_task_info = niftiinfo(fullfile(file_in.folder,file_in.name));
        main_task_data = niftiread(main_task_info);
        
        % trim scan, extract only main task phase
        main_start = behav_regressor_glm{1,4}{1,i}.scan_num(1);
        main_end = behav_regressor_glm{1,4}{1,i}.scan_num(2);
        extracted_data = main_task_data(:, :, :, main_start:main_end);
        
        % create nifti file _ main task only
        main_onset = main_task_info.PixelDimensions(4) * (main_start-1); % PixelDimensions(4) == TR을 의미
        main_task_info.ImageSize = size(extracted_data);
        current_nifti = fullfile(path_out_1st,[c_sbj '_main_task_trimmed.nii']);
        niftiwrite(extracted_data, current_nifti, main_task_info);
        
        
        %% 5) 1st GLM - Make multiple regressor
        % event regressors
        names = behav_regressor_glm{1, 4}{1, i}.regress_name;
        onsets = behav_regressor_glm{1, 4}{1, i}.regress_onset;
        durations = behav_regressor_glm{1, 4}{1, i}.regress_duration;
        save(fullfile(condition_path,[c_sbj, '_multiple_conditions.mat']), 'names', 'onsets', 'durations');
        
        % movement regressors
        movereg_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        conf_file = dir(fullfile(path_in_1st, c_sbj,'func', ['*task-OCAT', '*confounds*.tsv']));   % find the file
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
        
        %% 6) 1st GLM - fMRI model specification
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_spec.dir = {current_beta_out};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {current_nifti};
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
        
        %% 7) 1st GLM - Model estimation
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(current_beta_out, 'SPM.mat') };
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', matlabbatch)
        
        %% 8) 1st GLM - Contrast manager
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat  = {fullfile(current_beta_out, 'SPM.mat') };
  
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'correct';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0.5 0 0.5 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'incorrect';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 0.5 0 0.5];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Hit';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'correct_rejection';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Corr-Incorr';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0.5 -0.5 0.5 -0.5];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        if sbj_id_list(i)==7
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 0 1 0];
        end
        if sbj_id_list(i)==25
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 -1 1 0];
        end
        if sbj_id_list(i)==28
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 1 -1 0];
        end
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'Incorr-Corr';
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [-0.5 0.5 -0.5 0.5];
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'Hit-Corr_reject';
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [1 0 -1 0];
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'Corr_reject-Hit';
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [-1 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        spm_jobman('run', matlabbatch)
        
    end
end
%% %%%%%% 2nd level GLM! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_2nd_glm ==1
    clear matlabbatch;
    % make loop for each conditions
    batch_save_file=cell(numel(contrast{1,2}), 4);
    for ct_con=1:numel(contrast{1,2})
        % set path
        cond_path_out = fullfile(path_out_2nd,contrast{1,2}{ct_con}); mkdir(cond_path_out);
        % set files
        current_scans = char(arrayfun(@(i) [path_out_1st, sprintf('%s%02d','\sub-', i), '\',sprintf('%s%04d', 'con_',ct_con),'.nii,1'], sbj_id_list, 'UniformOutput', false));
        
        %% 9) 2nd GLM - Factorial design specification
        matlabbatch=[];
        % Directory
        matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(cond_path_out);
        % Design
        % Scans
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(string(current_scans));
        % Covariates
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        % Multiple covariates
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        % Masking
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        % Global calculation
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        % Global normalisation
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        % Normalisation
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        %% Model Estimation (2nd-level RFX)
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        %% Contrast Manager (2nd-level RFX)
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        % Name
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = contrast{1,2}{ct_con};
        % Weights vector
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        % Replicate over sessions (default; Don't replicate)
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        % Delete existing contrasts (Yes)
        matlabbatch{3}.spm.stats.con.delete = 1;
        
        %% Results Report
        if strcmp(threshold_type,'FWE')
            threshold = 0.05;
        elseif strcmp(threshold_type,'none')
            threshold = 0.001;
        elseif strcmp(threshold_type,'FDR')
            threshold = 0.05;
        end
        
        matlabbatch{4}.spm.stats.results.spmmat = {[cond_path_out '\SPM.mat']};
        matlabbatch{4}.spm.stats.results.conspec.titlestr = contrast{1,2}{ct_con};
        matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{4}.spm.stats.results.conspec.threshdesc = threshold_type;   % 'FWE' or 'FDR' or 'none'
        matlabbatch{4}.spm.stats.results.conspec.thresh = threshold;   % GUI로는 0.005로 했음
        matlabbatch{4}.spm.stats.results.conspec.extent = 0;    % Nah.m에서는 05로 했음
        matlabbatch{4}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
        matlabbatch{4}.spm.stats.results.units = 1;
        matlabbatch{4}.spm.stats.results.print = true;
        
        spm_jobman('run', matlabbatch);
        
        batch_save_file(ct_con,:)=matlabbatch(:);
    end
    
    %% Run & save batch
    batch_save_dir = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
    save([batch_save_dir '\RFX_final.mat'], 'batch_save_file');
    
end
