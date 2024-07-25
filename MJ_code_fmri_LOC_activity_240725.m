% 살제 activity 값들을 가져와서 mask==1인 부분만 mean을 낸 후
% 그 mean이 하나의 dot이 되어서 피험자들 모두를 bar plot 그리기
clear;clc;close
addpath('C:\Users\leelab\Documents\MATLAB\spm12');
%% roi mask 만드는 법 (using pickAtlas)
% wfu_pickatlas;

% set ROI mask
roi_masks = {'D:\sorin\data\rLOC_MNI.nii'};

mni=niftiread('C:\MRIcroGL\Resources\standard\rmni152.nii');

roi_mask.pfc=logical(niftiread(roi_masks{1}));
roi_mask.loc=logical(niftiread(roi_masks{2}));
roi_mask.hpc=logical(niftiread(roi_masks{3}));
%%
path_nii='D:\sorin\data\spm_ocat_bids\spmprep_final';
activity_scan={};all_values = [];
for sbj_i=1:numel(sbj_id_list_38)
    sbj_n = sbj_id_list_38(sbj_i);
    c_sbj=sprintf('sub_%d',sbj_n);
    func_nii=dir(fullfile(path_nii,c_sbj,'func','swraf*.nii'));
    activity_scan{1,sbj_i} = nan(1,height(func_nii));
    
    for shot=1:height(func_nii)
        curr_nii=niftiread(fullfile(func_nii(1).folder,func_nii(shot).name));
        %         current_pfc = curr_nii(pfc);
        current_loc = curr_nii(roi_mask.loc);
        activity_scan{1,sbj_i}(shot)=mean(current_loc,"all");
        all_values = [all_values, activity_scan{1, sbj_i}(shot)];
    end
end


%% normalization
mean_all = mean(all_values, 'omitnan');
std_all = std(all_values, 'omitnan');

for sbj_i=1:numel(sbj_id_list_38)
    
    activity_scan{1, sbj_i} = (activity_scan{1, sbj_i} - mean_all) / std_all;
end

%% 가장 기본적인 plot 형태:
figure()
plot(1:height(func_nii),activity_scan{1, sbj_i})
%%%%% 모든 피험자 동일한 설정으로 진행하기!!! %%%
% plot의 y축은 [-1.5:1.5]로 지정하고, 피험자마다 plot을 그려서 사진으로 저장해보세요
% figure size는 가로로 길게, 알아보기 쉽게 설정하기!(figure함수로 지정하기)
% plot의 선 굵기도 적당히 조절하기
% 검은화면 -> 흰화면 혹은 흰화면 -> 검은화면 바뀌는 scan_num을 찾아서 세로 선 그려주기 (x=## plot)






%% scan_number 참고!
load('regressors_GLM_0717.mat')
for sbj_i = 1:numel(sbj_id_list_38)
    sbj_n = sbj_id_list_38(sbj_i);
    c_sbj=sprintf('sub_%d'sbj_n);disp(c_sbj);
    
    main_scan_num{sbj_i}=reg_0717{1, sbj_i}.main.ver_1.scan_num;
    pre_ODT_scan_num{sbj_i}=reg_0717{1, sbj_i}.ODT.pre.scan_num;
    post_ODT_scan_num{sbj_i}=reg_0717{1, sbj_i}.ODT.post.scan_num;
end



