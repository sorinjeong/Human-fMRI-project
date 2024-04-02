% batch_newLSS

% India: Single trial regressor
% April, 2023

% General settings
% tStart = tic;
% load subject list
% load('D:\Delta_T_Analysis\Behavioral\Subject_data_final\Final_sbj_run_list.mat');
% sbj_run_list(2,:) = [];
% subjects = sbj_run_list;
% subject = {'709'};
% subject = '709';

% spmFolder = '/nfs/ep2/AX/first_levels/00_MONTH/';
% spmSubFolder = '/MTU/func_an_SPM8/';
%
% outFolder = '/home/tsalo/lssForUnivariate/';
% outSubFolder = '/';

% LSS settings
% includeConditions = {'Asso' 'TOJ_close' 'TOJ_far'};
% settings.model = 2;            % 1- Rissman, 2- LSS
% settings.useTempFS = 0;        % 0- do not use temporary files, 1- use temporary files
% settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% LSS settings
% includeConditions = {'CueA' 'CueB'};
% settings.model = 2;            % 1- Rissman, 2- LSS
% settings.useTempFS = 0;        % 0- do not use temporary files, 1- use temporary files
% settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings (not used in my analysis)
% rois = load('/nfs/cntracs/lssForKim/gt_rois_final.mat');
% settings.fConnType = 'roi2roi'; % seed2voxel or roi2roi
%
% for iSubj = 1:length(subjects)
%     spmDir = [spmFolder subjects{iSubj} spmSubFolder];
%     outDir = [outFolder subjects{iSubj} outSubFolder];
%
%     images = lssGenerateBetasSpm(subjects{iSubj}, spmDir, outDir, includeConditions, settings);
%     lssCorrelation(images, rois.rois, settings);
% end


%% Set path & execute lssGenerateBetasSPM function
% for extracting parameter estimates from each ROIs for each trials
% for extracting patterns from each ROIs for each trials


%% Run lssGenerateBetasSPM.m
clear; clc; close all;
%% Input!
main_or_ODT = 'main'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
threshold_type = 'FWE' ;% 'FWE' or 'FDR' or 'none'
date = '240330';
reg_mov_2or4 = 4; % if 4, moving regressor will contain [context(2)*order(first/remainder)], if 2, moving regressor will consider only context
phase = 'obj_show'; % main: 1) 'obj_show', 2) 'choice', 3) 'moving'
if strcmp(main_or_ODT,'ODT'); phase ='ODT';end

%set path
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');
root_path = '../data/data_fmri_glm';
addpath(genpath(root_path));

%% Get data
path_out_1st = fullfile(root_path, date,(main_or_ODT),(phase),'1st_Level');
load('regressors_GLM');
if reg_mov_2or4 == 4
    mov_ver = 'ver1';
    regressor = reg_for_glm_ver1;
elseif reg_mov_2or4 == 2
    mov_ver = 'ver2';
    regressor = reg_for_glm_ver2;
end
sbj_id_list
regressor

% preprocessed_path = 'C:\Users\zyeon\SynologyDrive\OCAT\data\brain_preprocessed';
% behavior_path = 'C:\Users\zyeon\SynologyDrive\OCAT\data\behavior';
% output_path = 'C:\Users\zyeon\SynologyDrive\OCAT\Results_new_smoothing\GLM\1stLevel';
%
% load(fullfile(behavior_path, 'regressor_all.mat'));
% sbj_id_list
% behav_regressor_glm
%%
% oc_list_name = {'OC_AssoObj','OC_NAObj', 'OC_correct','OC_incorrect','OC_AssoObj_Corr', 'OC_AssoObj_Incorr', 'OC_NAObj_Corr', 'OC_NAObj_Incorr', 'OC_AssoResp', 'OC_NAResp'};

% glm
% regress_num_standard = [8,8,10,8]; %%%%%%%%%%%

for sbj_i = 1:length(sbj_id_list)
    delete(gcp('nocreate'))

    % change this for different Design Matrix
    spmDir = fullfile(path_out_1st, sprintf('sub-%.2d',sbj_id_list(sbj_i)));
    outDir = fullfile(path_out_1st, 'single_trial', sprintf('sub-%.2d',sbj_id_list(sbj_i)));

    % spmDir = ['C:\Users\zyeon\SynologyDrive\OCAT\Results_new\GLM\1stLevel\glm' num2str(glm_i) '\SUB' num2str(sbj_id_list(sbj_i))];
    % outDir = ['C:\Users\zyeon\SynologyDrive\OCAT\Results_new\GLM\1stLevel\SinggleTrial\glm' num2str(glm_i) '\SUB' num2str(sbj_id_list(sbj_i))];

    % LSS settings
    % includeConditions = behav_regressor_glm{glm_i}{sbj_i}.regress_name;
    %         includeConditions = includeConditions(ismember(includeConditions, oc_list_name));
    includeConditions = regressor{1,sbj_i}.((main_or_ODT)).regress_name;


    settings.model = 2;            % 1- Rissman, 2- LSS
    settings.useTempFS = 0;        % 0- do not use temporary files, 1- use temporary files
    settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

    lssGenerateBetasSpm(sprintf('sub-%.2d',sbj_id_list(sbj_i)), spmDir, outDir, includeConditions, settings)

end  %-- end of for ct_sub




