% -- Modified by Sorin Jeong, February 2024

% Result Report
% Random Effect Analysis (2nd-level Model for Between-Subject Analysis
clear; clc; close all;
%set path
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');
root_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_glm';
% addpath(genpath(root_path));

%% Input!!
add_smoothing = 0;% if you want to add a smoothing step, put 1
run_1st_glm = 1; % if you want to run the 1st glm step, put 1
run_2nd_glm = 0; % if you want to run the 2nd glm step, put 1
main_or_ODT = 'ODT'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
threshold_type = 'FWE' ;% 'FWE' or 'FDR' or 'none'
date = '240512_new';

%% 1) Defining pathway and subjects
% load regressors
load('regressors_ODT_4obj_0513');
    regressor = odt_reg;
% sbj_id_list=sort([sbj_id_list 7]);

% set-up directory
path_in_1st = fullfile('../data','new_ocat_bids','derivatives'); % fmriprep's output file (before smoothing)

    condition_path = fullfile(root_path, date,(main_or_ODT),'data_for_glm');mkdir(condition_path);addpath(condition_path);
    path_out_1st = fullfile(root_path, date,(main_or_ODT),'1st_Level');
if ~exist(path_out_1st,"dir"); mkdir(path_out_1st);end

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


        %% 4) 1st GLM - Extract Scan
        % read 4D-scan
        file_in = dir(fullfile(sbj_dir,'func',[c_sbj '*preproc_bold.nii.gz']));
        nii_info = niftiinfo(fullfile(file_in.folder,file_in.name));
        nii_data = niftiread(nii_info);

        % trim scan, extract only main task phase

            pre_scan_start = regressor{1,i}.(main_or_ODT).pre_scan_num(1);
            pre_scan_end = regressor{1,i}.(main_or_ODT).pre_scan_num(2);
            post_scan_start = regressor{1,i}.(main_or_ODT).post_scan_num(1);
            post_scan_end = regressor{1,i}.(main_or_ODT).post_scan_num(2);
            scan_range = [pre_scan_start:pre_scan_end post_scan_start:post_scan_end];
        
        if scan_range(end)>nii_info.ImageSize(end)
            scan_range(end) = nii_info.ImageSize(end);
        end
        extracted_data = nii_data(:, :, :, scan_range);


        % create nifti file
        nii_info.ImageSize = size(extracted_data);
        current_nifti = fullfile(path_out_1st,[c_sbj '_' (main_or_ODT) '_trimmed.nii']);
        niftiwrite(extracted_data, current_nifti, nii_info);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 5) 1st GLM - Make multiple regressor
        % event regressors
        names = regressor{1,i}.((main_or_ODT)).regress_name;
        onsets = regressor{1,i}.((main_or_ODT)).regress_onset;
        durations = regressor{1,i}.((main_or_ODT)).regress_duration;
        save(fullfile(condition_path,[c_sbj, 'multi_regressor_ODT','.mat']), 'names', 'onsets', 'durations');


        % %% contrasts!!
        % contrast{1,1}='name: ';
        % contrast{2,1}='weight: ';
        % 
        % % ODT - contrast
        % cont_name_typ3 = {'pre-ODT'; 'post-ODT'; 'post-pre'; 'pre-forest'; 'pre-city'; 'post-forest'; 'post-city'; 'forest_post-pre'; 'city_post-pre'; 'target'; 'target_post-pre'};
        % cont_weight_typ3={[0.25 0.25 0.25 0.25]; [0 0 0 0 0.25 0.25 0.25 0.25]; [-0.25 -0.25 -0.25 -0.25 0.25 0.25 0.25 0.25]; [0.5 0.5 0 0]; [0 0 0.5 0.5]; ...
        %     [0 0 0 0 0.5 0.5]; [0 0 0 0 0 0 0.5 0.5]; [-0.5 -0.5 0 0 0.5 0.5]; [0 0 -0.5 -0.5 0 0 0.5 0.5]; [0 0 0 0 0 0 0 0 0.5 0.5]; [0 0 0 0 0 0 0 0 -0.5 0.5]};

        % %% phase of Interest _ Contrast
        % 
        %     contrast{1,2}=cont_name_typ3;contrast{2,2}=cont_weight_typ3;
        %     if isnan(names{end})
        %         names(end) = []; onsets(end) = []; durations(end) = [];
        %     end
        %     durations{11, 1} = 4;
        % 
        % 
        % save(char(fullfile(root_path,date,(main_or_ODT),['contrast_ODT','.mat'])),"contrast");


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

        save(fullfile(condition_path, [c_sbj, '_movements.mat']), 'R', 'R_names');

        %% 6) 1st GLM - fMRI model specification
        addpath('C:\Users\Leelab\Documents\MATLAB\spm12')

        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_spec.dir = {current_beta_out};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {current_nifti};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(condition_path,[c_sbj, strcat('multi_regressor_ODT','.mat')])};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(condition_path, [c_sbj, '_movements.mat'])};
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

        % %% 8) 1st GLM - Contrast manager
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
end
