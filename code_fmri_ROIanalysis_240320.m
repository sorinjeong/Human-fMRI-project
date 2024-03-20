% 지연이가 알려준대로 1st level에서 개별적으로 뽑은 beta값들을 가져와서 mask==1인 부분만 mean을 낸 후
% 그 mean이 하나의 dot이 되어서 피험자들 모두를 bar plot 그리기
clear;clc;close
%% roi mask 만드는 법 (using pickAtlas)

% addpath(genpath('C:\Users\leelab\Documents\MATLAB\spm12\toolbox\WFU_PickAtlas_3.0.5b'));
% wfu_pickatlas;

%% set ROI mask
% mni_hpc_mask_file = 'D:\fMRI\OCAT_DIR\data\atlas\Standard-space\mni152_hippocampal_area.nii';
% mni_hpc_coord=mask2coord(mni_hpc_mask_file);

which_roi = 1:9; %put roi number you want; *FYI about roi number: see roi_masks variable
roi_masks = {'L_HPC';'L_PHG';'L_hipp_area';'R_HPC';'R_PHG';'R_hipp_area';'bi_HPC';'bi_PHG';'bi_hipp_area'};

% roi_file = {'mni152_hippocampal_area', 'mni152_R_HPC', 'mni152_L_HPC','mni152_R_PHG','mni152_L_PHG','mni152_ant_HPC','mni152_post_HPC'};
roi_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\atlas\Standard-space\hpc';
roi_files = dir(fullfile(roi_path,'mni152_*.nii'));

%% Input!!!
main_or_ODT = 'ODT'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
% phase of Interest
% main: 1) 'obj_show', 2) 'choice', 3) 'obj_ITI', 4) 'moving'
% ODT: 'ODT'
phase = 'ODT'; % input
if strcmp(main_or_ODT,'ODT'); phase ='ODT';end

%% set path and load variables
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');

path_root = '../data/240320';
path_out_1st = fullfile(path_root,(main_or_ODT),(phase),'1st_Level');
path_out_2nd = fullfile(path_root,(main_or_ODT),(phase),'2nd_Level');

load(fullfile(path_root,(main_or_ODT),(phase), 'contrast_v1.mat'));
load('reg_for_glm_mov');

%% define the contrast names
cont_name_typ1 = {'correct','incorrect','Hit', 'miss','flase_positive','correct_rejection','Corr-Incorr','Incorr-Corr','Hit-Corr_reject','Corr_reject-Hit','Hit-miss','Corr_reject-false_pos'};
cont_name_typ2 = {'forest', 'city', 'forest-city', 'first_seen', 'navigation','first_seen_vs_navigate'};
cont_name_typ3 = {'pre-ODT', 'post-ODT', 'post-pre', 'pre-forest', 'pre-city', 'post-forest', 'post-city', 'forest_post-pre', 'city_post-pre', 'target', 'target_post-pre'};

% Determine the contrast based on main_or_ODT and phase
if strcmp(main_or_ODT, 'main') && ~strcmp(phase, 'moving')
    con_name = cont_name_typ1;
elseif strcmp(main_or_ODT, 'main') && strcmp(phase, 'moving')
    con_name = cont_name_typ2;
elseif strcmp(main_or_ODT, 'ODT')
    con_name = cont_name_typ3;
end

%% roi 별 mean값 구하기
conditions_roi={};
for ct_con=1:numel(con_name)
    disp(con_name{ct_con})
    current_con={}; % 여기서 초기화하기!!!
    path_condition = fullfile(path_root,'ROI_240320',(main_or_ODT),(phase));mkdir(path_condition);

    T = array2table(nan(numel(sbj_id_list), numel(roi_masks)), 'VariableNames', roi_masks);
    for r=which_roi
        current_roi = fullfile(roi_path,roi_files(r).name);
        roi_coord = mask2coord(current_roi);
        for n_sbj=1:numel(sbj_id_list)
            c_sbj = sprintf('sub-%.2d', sbj_id_list(n_sbj));
            
            %% get the list of contrast's beta
            sbj_beta = dir(fullfile(path_out_1st, c_sbj, ['con_*_' con_name{ct_con} '.nii']));
            % sbj_beta = fullfile(path_out_1st, c_sbj, sprintf('%s%04d%s', 'con_',ct_con,'.nii'));
            if ~isempty(sbj_beta)
                roi_beta = niftiread(fullfile(sbj_beta(1).folder,sbj_beta(1).name));
                A=[];
                for i=1:size(roi_coord,1)
                    x=roi_coord(i,1); y=roi_coord(i,2); z=roi_coord(i,3);
                    A(i)=roi_beta(x,y,z);
                    roi_mean(n_sbj,1) = nanmean(A);
                end
            end
        end
        current_con{r,1} = con_name{ct_con};
        current_con{r,2} = roi_masks{r};
        current_con{r,3} = roi_mean;
        T.(roi_masks{r}) = roi_mean;
    end
    %     n=numel(which_roi)*(ct_con-1)+1;
    %     conditions_roi{n:n+numel(which_roi),:} = current_con;
    conditions_roi = vertcat(conditions_roi, current_con);
    writetable(T, fullfile(path_condition,'conditions_roi_table.xlsx'), 'Sheet', con_name{ct_con});
    save(fullfile(path_condition,strcat(con_name{ct_con},'.mat')),"current_con");
end

save(fullfile(path_condition,'conditions_roi.mat'),"conditions_roi");

