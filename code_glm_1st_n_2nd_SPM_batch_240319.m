% -- Modified by Sorin Jeong, February 2024

% Result Report
% Random Effect Analysis (2nd-level Model for Between-Subject Analysis
clear; clc; close all;
%set path
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');
root_path = '../data';
addpath(genpath(root_path));

%% Input!!
add_smoothing = 0;% if you want to add a smoothing step, put 1
run_1st_glm = 1; % if you want to run the 1st glm step, put 1
run_2nd_glm = 0; % if you want to run the 2nd glm step, put 1
main_or_ODT = 'main'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
threshold_type = 'FWE' ;% 'FWE' or 'FDR' or 'none'
date = '240319_new3';
% isROI = 1; % if you want to run the ROI analysis, put 1, or if not, put 0 for whole brain analysis

%% phase of Interest
% main: 1) 'obj_show', 2) 'choice', 3) 'obj_ITI', 4) 'moving'
% ODT: 'ODT'

phase = 'obj_show'; % input


%% 1) Defining pathway and subjects

% set-up directory
path_in_1st = fullfile(root_path,'data_fmri_bids','derivatives'); % fmriprep's output file (before smoothing)
path_out_1st = fullfile(root_path, date,(main_or_ODT),(phase),'1st_Level'); % save betas this folder
path_out_2nd = fullfile(root_path,date,(main_or_ODT),(phase),'2nd_Level'); % save 2nd lev glm results this folder
if ~exist(path_out_1st); mkdir(path_out_1st);end
if ~exist(path_out_2nd); mkdir(path_out_2nd);end

condition_path = fullfile(root_path, date,(main_or_ODT),(phase),'data_for_glm');mkdir(condition_path);addpath(condition_path);

% load regressors
load('reg_for_glm_mov');

% create subjects beta and contrast result
beta_contrast_table = cell(length(sbj_id_list), 3);

%% 2) Initialize SPM
spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;

%% 3) preprocessing - Spatial Smoothing
if run_1st_glm == 1
    for i = 1:length(sbj_id_list)
        c_sbj = sprintf('sub-%.2d', sbj_id_list(i)); disp(c_sbj)

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
        nii_info = niftiinfo(fullfile(file_in.folder,file_in.name));
        nii_data = niftiread(nii_info);

        % trim scan, extract only main task phase
        if strcmp((main_or_ODT),'main')
            scan_start = reg_for_glm{1,i}.(main_or_ODT).scan_num(1);
            scan_end = reg_for_glm{1,i}.(main_or_ODT).scan_num(2);
            scan_range = [scan_start:scan_end];
        elseif strcmp((main_or_ODT),'ODT')
            pre_scan_start = reg_for_glm{1,i}.(main_or_ODT).pre_scan_num(1);
            pre_scan_end = reg_for_glm{1,i}.(main_or_ODT).pre_scan_num(2);
            post_scan_start = reg_for_glm{1,i}.(main_or_ODT).post_scan_num(1);
            post_scan_end = reg_for_glm{1,i}.(main_or_ODT).post_scan_num(2);
            scan_range = [pre_scan_start:pre_scan_end post_scan_start:post_scan_end];
        end
        extracted_data = nii_data(:, :, :, scan_range);


        % create nifti file
        nii_info.ImageSize = size(extracted_data);
        current_nifti = fullfile(path_out_1st,[c_sbj '_' (main_or_ODT),(phase) '_trimmed.nii']);
        niftiwrite(extracted_data, current_nifti, nii_info);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 5) 1st GLM - Make multiple regressor
        % event regressors
        names = reshape(reg_for_glm{1,i}.((main_or_ODT)).regress_name',[],1);
        onsets = reshape(reg_for_glm{1,i}.(main_or_ODT).regress_onset',[],1);
        durations = reshape(reg_for_glm{1,i}.(main_or_ODT).regress_duration',[],1);

        % remove empty condition
        emptyIdx = cellfun(@(x) any(isempty(x)), onsets);

        % re-enter the contrast after removing conditions
        if sum(emptyIdx)>0
            disp('The following names had empty onsets:');
            disp(names(emptyIdx));

            emptySuffixes = cellfun(@(x) x(find(x=='_', 1, 'last')+1:end), names(emptyIdx), 'UniformOutput', false);
            names(emptyIdx) = []; onsets(emptyIdx) = []; durations(emptyIdx) = [];
            zero_cont = zeros(1,5-(length(emptyIdx)/3));
        else
            zero_cont = zeros(1,5);
        end

        %% contrasts!!
        contrast{1,1}='name: ';
        contrast{2,1}='weight: ';

        % main - contrast
        cont_name_typ1 = {'correct';'incorrect';'Hit'; 'miss';'flase_positive';'correct_rejection';'Corr-Incorr';'Incorr-Corr';'Hit-Corr_reject';'Corr_reject-Hit';'Hit-miss';'Corr_reject-false_pos'};
        cont_weight_typ1 = {[0.5 0 0 0.5]; [0 0.5 0.5 0]; [1 0 0 0]; [0 1 0 0]; [0 0 1 0]; [0 0 0 1]; [0.5 -0.5 -0.5 0.5]; [-0.5 0.5 0.5 -0.5]; [1 0 0 -1]; [-1 0 0 1]; [1 -1 0 0]; [0 0 -1 1]};
        if sum(emptyIdx)>0
            if contains(emptySuffixes,'hit')
                cont_name_typ1 = {'correct';'incorrect'; 'miss';'flase_positive';'correct_rejection';'Corr-Incorr';'Incorr-Corr';'Corr_reject-false_pos'};
                cont_weight_typ1 = {[0 0 0 1]; [0 0.5 0.5 0];  [0 1 0 0]; [0 0 1 0]; [0 0 0 1]; [0 -0.5 -0.5 1]; [0 0.5 0.5 -1]; [0 0 -1 1]};
            elseif contains(emptySuffixes,'miss')
                cont_name_typ1 = {'correct';'incorrect';'Hit';'flase_positive';'correct_rejection';'Corr-Incorr';'Incorr-Corr';'Hit-Corr_reject';'Corr_reject-Hit';'Corr_reject-false_pos'};
                cont_weight_typ1 = {[0.5 0 0 0.5]; [0 0 1 0]; [1 0 0 0]; [0 0 1 0]; [0 0 0 1]; [0.5 0 -1 0.5]; [-0.5 0 1 -0.5]; [1 0 0 -1]; [-1 0 0 1]; [0 0 -1 1]};

            elseif contains(emptySuffixes,'false')
                cont_name_typ1 = {'correct';'incorrect';'Hit'; 'miss'; 'correct_rejection';'Corr-Incorr';'Incorr-Corr';'Hit-Corr_reject';'Corr_reject-Hit';'Hit-miss'};
                cont_weight_typ1 = {[0.5 0 0 0.5]; [0 1 0 0]; [1 0 0 0]; [0 1 0 0];  [0 0 0 1]; [0.5 -1 0 0.5]; [-0.5 1 0 -0.5]; [1 0 0 -1]; [-1 0 0 1]; [1 -1 0 0]};

            elseif contains(emptySuffixes,'corr_rej')
                cont_name_typ1 = {'correct';'incorrect';'Hit'; 'miss';'flase_positive';'Corr-Incorr';'Incorr-Corr';'Hit-miss'};
                cont_weight_typ1 = {[1 0 0 0]; [0 0.5 0.5 0]; [1 0 0 0]; [0 1 0 0]; [0 0 1 0];  [1 -0.5 -0.5 0]; [-1 0.5 0.5 0]; [1 -1 0 0]};

            end
        end
        % main-type2: moving
        cont_name_typ2 = {'forest'; 'city'; 'forest-city'; 'first_seen'; 'navigation';'first_seen_vs_navigate'};
        cont_weight_typ2={[0.5 0 0.5 0]; [0 0.5 0 0.5]; [0.5 -0.5 0.5 -0.5]; [0.5 0.5 0 0]; [0 0 0.5 0.5]; [0.5 0.5 -0.5 -0.5]};

        % ODT - contrast
        cont_name_typ3 = {'pre-ODT'; 'post-ODT'; 'post-pre'; 'pre-forest'; 'pre-city'; 'post-forest'; 'post-city'; 'forest_post-pre'; 'city_post-pre'; 'target'; 'target_post-pre'};
        cont_weight_typ3={[0.25 0.25 0.25 0.25]; [0 0 0 0 0.25 0.25 0.25 0.25]; [-0.25 -0.25 -0.25 -0.25 0.25 0.25 0.25 0.25]; [0.5 0.5 0 0]; [0 0 0.5 0.5]; ...
            [0 0 0 0 0.5 0.5]; [0 0 0 0 0 0 0.5 0.5]; [-0.5 -0.5 0 0 0.5 0.5]; [0 0 -0.5 -0.5 0 0 0.5 0.5]; [0 0 0 0 0 0 0 0 0.5 0.5]; [0 0 0 0 0 0 0 0 -0.5 0.5]};

        %% phase of Interest _ Contrast

        % main
        if strcmp(phase,'obj_show')
            phase_idx = {1,1:4}; contrast{1,2}=cont_name_typ1; contrast{2,2}=cont_weight_typ1;
        elseif strcmp(phase,'choice')
            phase_idx = {2,1:4}; contrast{1,2}=cont_name_typ1;
            contrast{2,2}= cellfun(@(x) [ zero_cont x], cont_weight_typ1, 'UniformOutput', false);
        elseif strcmp(phase,'obj_ITI')
            phase_idx = {3,1:4}; contrast{1,2}=cont_name_typ1;
            contrast{2,2}= cellfun(@(x) [ zero_cont zero_cont x], cont_weight_typ1, 'UniformOutput', false);
        elseif strcmp(phase,'moving')
            phase_idx = {4,1:4};contrast{1,2}=cont_name_typ2;
            contrast{2,2}= cellfun(@(x) [ zero_cont zero_cont zero_cont x], cont_weight_typ2, 'UniformOutput', false);

            % ODT
        elseif strcmp(phase,'ODT')
            contrast{1,2}=cont_name_typ3;contrast{2,2}=cont_weight_typ3;
        end

        save(char(fullfile(root_path,date,(main_or_ODT),(phase),'contrast_v1.mat')),"contrast");

        save(fullfile(condition_path,[c_sbj, '_multiple_conditions.mat']), 'names', 'onsets', 'durations');

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
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; %autoregrssive ; AR(1) is the default option (by AHN)
        spm_jobman('run', matlabbatch)

        %% 7) 1st GLM - Model estimation
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(current_beta_out, 'SPM.mat') };
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', matlabbatch)

        %% 8) 1st GLM - Contrast manager
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat  = {fullfile(current_beta_out, 'SPM.mat') };
        for c = 1:length(contrast{1,2})
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = contrast{1,2}{c};
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = contrast{2,2}{c,:};
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';

        end
        matlabbatch{1}.spm.stats.con.delete = 1;

        spm_jobman('run', matlabbatch)


        %% rename the beta, contrast files for each subject
        beta_files = dir(fullfile(current_beta_out, 'beta_*.nii'));
        for b = 1:length(beta_files)
            old_name = fullfile(current_beta_out,beta_files(b).name);
            new_name = strcat(strrep(old_name, '.nii', '_'), names{b}, '.nii');
            copyfile(old_name, new_name);
        end

        con_files = dir(fullfile(current_beta_out, 'con_*.nii'));
        for cn = 1:length(con_files)
            old_name = fullfile(current_beta_out,con_files(cn).name);
            new_name = strcat(strrep(old_name, '.nii', '_'), contrast{1,2}{cn}, '.nii');
            copyfile(old_name, new_name);
        end
        beta_contrast_table{i,1} = c_sbj;
        beta_contrast_table{i,2} = dir(fullfile(current_beta_out, 'beta_*_*.nii'));
        beta_contrast_table{i,3} = dir(fullfile(current_beta_out, 'con_*_*.nii'));

    end
end
result_T = cell2table(beta_contrast_table, 'VariableNames', {'subject', 'condition','contrast'});
writetable(result_T,[path_out_1st '\result_T.xlsx']);

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
    save([batch_save_dir '\RFX_240314.mat'], 'batch_save_file');

end
