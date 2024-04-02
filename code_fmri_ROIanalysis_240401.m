% 지연이가 알려준대로 1st level에서 개별적으로 뽑은 beta값들을 가져와서 mask==1인 부분만 mean을 낸 후
% 그 mean이 하나의 dot이 되어서 피험자들 모두를 bar plot 그리기
clear;clc;close
%% roi mask 만드는 법 (using pickAtlas)

% addpath(genpath('C:\Users\leelab\Documents\MATLAB\spm12\toolbox\WFU_PickAtlas_3.0.5b'));
% wfu_pickatlas;

%% Input!!!
main_or_ODT = 'main'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
% phase of Interest
% main: 1) 'obj_show', 2) 'choice', 3) 'moving'
% ODT: 'ODT'
phase = 'moving'; % input
if strcmp(main_or_ODT,'ODT'); phase ='ODT';end
reg_mov_2or4 = 4; % if 4, moving regressor will contain [context(2)*order(first/remainder)], if 2, moving regressor will consider only context

%% set path and load variables
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');

path_root = '../data/data_fmri_glm/240330';
path_out_1st = fullfile(path_root,(main_or_ODT),(phase),'1st_Level');
path_out_2nd = fullfile(path_root,(main_or_ODT),(phase),'2nd_Level');
roi_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\atlas\local';addpath(roi_path);

load(fullfile(path_root,(main_or_ODT),(phase), 'contrast_ver1.mat'));
load('regressors_GLM');

%% define the contrast names
cont_name_typ1 = {'correct';'incorrect';'Corr-Incorr';'Incorr-Corr'};
cont_name_typ2_1 = {'forest'; 'city'; 'forest-city'; 'first_seen'; 'navigation';'first_seen_vs_navigate'};
cont_name_typ2_2 = {'forest'; 'city'; 'forest-city'};

cont_name_typ3 = {'pre-ODT', 'post-ODT', 'post-pre', 'pre-forest', 'pre-city', 'post-forest', 'post-city', 'forest_post-pre', 'city_post-pre', 'target', 'target_post-pre'};

% Determine the contrast based on main_or_ODT and phase
if strcmp(main_or_ODT, 'main') && ~strcmp(phase, 'moving')
    con_name = cont_name_typ1;
elseif strcmp(main_or_ODT, 'main') && strcmp(phase, 'moving') && reg_mov_2or4 == 4
    con_name = cont_name_typ2_1;
elseif strcmp(main_or_ODT, 'main') && strcmp(phase, 'moving') && reg_mov_2or4 == 2
    con_name = cont_name_typ2_2;
elseif strcmp(main_or_ODT, 'ODT')
    con_name = cont_name_typ3;
end

%% roi 별 mean값 구하기
conditions_roi={};roi_mean=[];
for ct_con=1:numel(con_name)
    disp(con_name{ct_con})
    current_con={}; % 여기서 초기화하기!!!
    path_condition = fullfile(path_root,'ROI_result',(main_or_ODT),(phase));mkdir(path_condition);
    %% set ROI mask

    for n_sbj=1:numel(sbj_id_list)
        load(fullfile(roi_path,strcat(string(n_sbj), '.mat')));
        current_mask = seg.seg_fit.hpc;
        % roi_masks = current_mask.roi_name_list;
        roi_masks = cellfun(@(x) strrep(x,'.','_'), current_mask.roi_name_list, UniformOutput=false);
        roi_files = current_mask.roi_list_nn;
        c_sbj = sprintf('sub-%.2d', sbj_id_list(n_sbj));
        T = array2table(nan(numel(sbj_id_list), numel(roi_masks)), 'VariableNames', roi_masks);
        for r=1:numel(roi_masks)
            current_roi = roi_files{r};
            %% get the list of contrast's beta
            sbj_beta = dir(fullfile(path_out_1st, c_sbj, strcat('con_*', string(ct_con), '.nii')));
            if ~isempty(sbj_beta)
                roi_beta = niftiread(fullfile(sbj_beta(1).folder,sbj_beta(1).name));
                ROI=roi_beta(current_roi);
                roi_mean = nanmean(ROI);
            end

            current_con{r,1} = con_name{ct_con};
            current_con{r,2} = roi_masks{r};
            current_con{r,3} = roi_mean;
            T.(roi_masks{r})(n_sbj) = roi_mean;
        end
    end
    conditions_roi = vertcat(conditions_roi, current_con);
    writetable(T, fullfile(path_condition,'conditions_roi_table.xlsx'), 'Sheet', con_name{ct_con});
    save(fullfile(path_condition,strcat(con_name{ct_con},'.mat')),"current_con");
end

save(fullfile(path_condition,'conditions_roi.mat'),"conditions_roi");

