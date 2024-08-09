% 살제 activity 값들을 가져와서 mask==1인 부분만 mean을 낸 후
% 그 mean이 하나의 dot이 되어서 피험자들 모두를 plot 그리기
% clear;clc;close
addpath('C:\Users\Leelab\Documents\MATLAB\spm12');

cd('D:\MJ\MJ_share')
% scan_num 열 추가
sheet_names=sheetnames("event_TR.xlsx");
numsheets=numel(sheet_names);
for i= 1: numsheets
    c_sbj=sprintf('sub-%02d',i);
    curr_T= readtable("event_TR.xlsx","Sheet",c_sbj);
    %%
    first_scan=find(curr_T.Var2==0,1,'last');
    curr_T.scan_num=[zeros(first_scan-1,1);(1:height(curr_T)-first_scan+1)'];
    writetable(curr_T,"event_TR.xlsx","Sheet",c_sbj)
end
%% roi mask 만드는 법 (using pickAtlas)
% wfu_pickatlas;

% set ROI mask
roi_masks = {'rLOC_MNI.nii'};
mni=niftiread('rmni152.nii');
roi_mask.loc=logical(niftiread(roi_masks{1}));

%%
load('regressors_GLM_0717.mat')
path_nii='D:\MJ\spmprep_final';
activity_scan={};all_values = [];
for sbj_i=1:numel(sbj_id_list_38)
    sbj_n = sbj_id_list_38(sbj_i);
    c_sbj=sprintf('sub_%d',sbj_n);disp(c_sbj)
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

    activity_scan_norm{1, sbj_i} = (activity_scan{1, sbj_i} - mean_all) / std_all;
end

%% 가장 기본적인 plot 형태:

% figure('Position',get(0, 'Screensize'))
for sbj_i=1:2%numel(sbj_id_list_38)
    % for sbj_i=[3,5,14]
    figure('Position',[300,300,1300,500])
    sbj_n = sbj_id_list_38(sbj_i);
    c_sbj=sprintf('sub_%d',sbj_n);disp(c_sbj);


    main_scan_num{sbj_i}=reg_0717{1, sbj_i}.main.scan_num;
    pre_ODT_scan_num{sbj_i}=reg_0717{1, sbj_i}.ODT.pre.scan_num;
    post_ODT_scan_num{sbj_i}=reg_0717{1, sbj_i}.ODT.post.scan_num;

    curr_pre=pre_ODT_scan_num{sbj_i};
    curr_main=main_scan_num{sbj_i};
    curr_post=post_ODT_scan_num{sbj_i};

    func_nii=dir(fullfile(path_nii,sprintf('sub_%d',sbj_id_list_38(sbj_i)),'func','swraf*.nii'));

    %plotting
    % subplot(5,7,sbj_i) %--> ylim 설정할때는 이거 키고, 아닐때는 %로 각주달기! / 위치는 항상 plot 위에 두어야 함!

    plot(1:height(func_nii),activity_scan{1, sbj_i},'Linewidth',3)
    %% figure 창 이름 설정하기!
    hold on
    title(['\fontsize{16}',sprintf('sub-%d',sbj_id_list_38(sbj_i)),' activity LOC'])
    box off
    % ylim([-2.5,3])

    %pre
    pre_main_x=[curr_pre(end):curr_main(1)];
    %all_subs
    % pre_main_y=repmat(1.8,size(pre_main_x));
    %individual
    pre_main_y=repmat((max(activity_scan{1, sbj_i})+0.05),size(pre_main_x));

    plot(pre_main_x,pre_main_y)
    line(pre_main_x,pre_main_y,'Linewidth',2)
    xline([curr_pre(end)+3,curr_main(1)+3])

    %post
    main_post_x=[curr_main(end):curr_post(1)];
    %all_subs
    % main_post_y=repmat(1.8,size(main_post_x));
    %individual
    main_post_y=repmat((max(activity_scan{1, sbj_i})+0.05),size(main_post_x));

    plot(main_post_x,main_post_y)
    line(main_post_x,main_post_y,'Linewidth',2)
    xline([curr_main(end)+3,curr_post(1)+3])

    % exportgraphics(gcf,([sprintf('sub-%d',sbj_id_list_38(sbj_i)),' activity LOC','.png']))
end



%%%%% 모든 피험자 동일한 설정으로 진행하기!!! %%%
% plot의 y축은 [-1.5:1.6]로 지정하고, 피험자마다 plot을 그려서 사진으로 저장해보세요
% figure size는 가로로 길게, 알아보기 쉽게 설정하기!(figure함수로 지정하기)
% plot의 선 굵기도 적당히 조절하기
% 검은화면 -> 흰화면 혹은 흰화면 -> 검은화면 바뀌는 scan_num을 찾아서 세로 선 그려주기 (x=## plot)


% y=1.5 지점에


%% scan_number 참고!




