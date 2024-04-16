

%% set parameters
clear;clc;cd('B:\OCAT_DIR\code\')
load('regressors_GLM_0405')
beta_root = 'B:\OCAT_DIR\data\data_fmri_glm\240403\main\obj_show\1st_Level\single_trial';

%% input!!
cond_or_trial = 'cond'; %cond별로 뽑을래, trial별로 뽑을래?; '1. 'cond', 2. 'trial'



%% condition별 parsing
cond = {'hit', 'miss', 'false', 'corr_rej'};
% cond = {'hit'};
% folder name
cond_beta=struct;dir_info=struct;
for i_sbj = 1: numel(sbj_id_list)
    n_sbj = sbj_id_list(i_sbj);
    c_sbj = sprintf('sub-%.2d',n_sbj);
    for c_i = 1: numel(cond)
        dir_info = dir(fullfile(beta_root,c_sbj,'betas','naming',strcat('*',(cond{c_i}),'*'),'beta_0001.nii'));
        cond_beta.(cond{c_i}){i_sbj,1} = arrayfun(@(x)  fullfile(dir_info(x).folder, dir_info(x).name), 1:height(dir_info), 'UniformOutput', false);
    end
end

%% Trial별 parsing
% corr 중에서 trial 줄세우기
trial = {'corr','incorr','total'};
% trial = {'total'};
trial_beta=struct;
for i_sbj = 1: numel(sbj_id_list)
    n_sbj = sbj_id_list(i_sbj);
    c_sbj = sprintf('sub-%.2d',n_sbj);
    for t_i = 1: numel(trial)
        dir_info = dir(fullfile(beta_root,c_sbj,'betas','naming',strcat('*_',(trial{t_i}),'_*')));
        if strcmp(trial{t_i},'total'); dir_info = dir(fullfile(beta_root,c_sbj,'betas','naming','*corr*'));end
        trial_nums = cellfun(@(x) str2double(regexp(x, '\d+$', 'match')), {dir_info.name});
        % trial 번호로 정렬
        [~, sort_idx] = sort(trial_nums);
        sort_names = {dir_info.name}';sort_names = sort_names(sort_idx);
        sort_folders = {dir_info.folder}';sort_folders=sort_folders(sort_idx);
        trial_beta.(trial{t_i}){i_sbj,1} = cellfun(@(x, y) fullfile(y, x, 'beta_0001.nii'), sort_names, sort_folders, 'UniformOutput', false);
        trial_beta.(trial{t_i}){i_sbj,2} = sort_idx';

    end
end

%% %%%%%%%%%%%%%  ROI analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_roi = [1 2 4 5 7 8]; %put roi number you want; *FYI about roi number: see roi_masks variable
roi_masks = {'L_HPC';'L_PHG';'L_hipp_area';'R_HPC';'R_PHG';'R_hipp_area';'bi_HPC';'bi_PHG';'bi_hipp_area'};

roi_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\atlas\Standard-space\hpc';
roi_files = dir(fullfile(roi_path,'mni152_*.nii'));

%% ROI_Input!!!
main_or_ODT = 'main'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
% phase of Interest
% main: 1) 'obj_show', 2) 'choice', 3) 'moving'
% ODT: 'ODT'
phase = 'obj_show'; % input
if strcmp(main_or_ODT,'ODT'); phase ='ODT';end
reg_mov_2or4 = 4; % if 4, moving regressor will contain [context(2)*order(first/remainder)], if 2, moving regressor will consider only context
if reg_mov_2or4 == 4; mov_ver = 'ver1';
elseif reg_mov_2or4 == 2; mov_ver = 'ver2';
end

%% set path and load variables

path_root = '../data/data_fmri_glm/240403';
if strcmp(phase,'moving')
    path_out_1st = fullfile(path_root,(main_or_ODT),(phase),(mov_ver),'1st_Level');
else
    path_out_1st = fullfile(path_root,(main_or_ODT),(phase),'1st_Level');
end

load(fullfile(path_root,(main_or_ODT),(phase), ['contrast_',mov_ver,'.mat']));


%% roi 별 mean값 구하기

%% cond
if strcmp(cond_or_trial,'cond')
    con_name=cond;
    beta=cond_beta;

    conditions_roi={};
    for con_i=1:numel(con_name)
        disp([cond_or_trial ' : ' con_name{con_i}])
        current_con={}; % 여기서 초기화하기!!!
        path_condition = fullfile(path_root,'ROI_result',(main_or_ODT),(phase),(mov_ver),(cond_or_trial));
        mkdir(path_condition);

        T = array2table(nan(numel(sbj_id_list), numel(which_roi)), 'VariableNames', roi_masks(which_roi));
        roi_mean=[];
        for r=1:numel(which_roi)
            current_roi = fullfile(roi_path,roi_files(which_roi(r)).name);
            roi_coord = mask2coord(current_roi);
            for i_sbj = 1: numel(sbj_id_list)
                sbj_beta = beta.(con_name{con_i}){i_sbj,1};
                beta_roi_mean=[];
                for b=1:numel(sbj_beta)
                    roi_beta = niftiread(sbj_beta{b});
                    A=[];
                    for i=1:size(roi_coord,1)
                        x=roi_coord(i,1); y=roi_coord(i,2); z=roi_coord(i,3);
                        A(i)=roi_beta(x,y,z);
                        beta_roi_mean(b,1) = nanmean(A);
                    end

                    roi_mean(i_sbj,1) = nanmean(beta_roi_mean);
                end
            end
            current_con{r,1} = con_name{con_i};
            current_con{r,2} = roi_masks{which_roi(r)};
            current_con{r,3} = roi_mean;
            T.(roi_masks{which_roi(r)}) = roi_mean;
        end
        conditions_roi = vertcat(conditions_roi, current_con);
        writetable(T, fullfile(path_condition,strcat(phase,'_roi_table.xlsx')), 'Sheet', con_name{con_i});
        save(fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.mat')),"current_con");
    end

    save(fullfile(path_condition,strcat(phase,'_conditions_roi.mat')),"conditions_roi");

end



%% trial ; corr and incorr

if strcmp(cond_or_trial,'trial')
    con_name=trial(1:2);
    beta=trial_beta;
    trials=1:32;


    conditions_roi={};

    for con_i=1:numel(con_name)
        disp([cond_or_trial ' : ' con_name{con_i}])
        current_con={}; % 여기서 초기화하기!!!
        path_condition = fullfile(path_root,'ROI_result',(main_or_ODT),(phase),(mov_ver),(cond_or_trial));
        mkdir(path_condition);

        for t=trials
            current_trial={}; % 여기서 초기화하기!!!

T = array2table(nan(numel(sbj_id_list), numel(which_roi)), 'VariableNames', roi_masks(which_roi));
            roi_mean=[];
            trial_roi_mean=zeros(1,32);

            for r=1:numel(which_roi)
                current_roi = fullfile(roi_path,roi_files(which_roi(r)).name);
                roi_coord = mask2coord(current_roi);
                for i_sbj = 1: numel(sbj_id_list)
                    sbj_beta = beta.(con_name{con_i}){i_sbj,1};
                    trial_idx = beta.(con_name{con_i}){i_sbj,2};
                    if strcmp(con_name{con_i},'total')
                        trial_idx=trials; end

                    for b=1:numel(sbj_beta)
                        roi_beta = niftiread(sbj_beta{b});

                        A=[];
                        for i=1:size(roi_coord,1)
                            x=roi_coord(i,1); y=roi_coord(i,2); z=roi_coord(i,3);
                            A(i)=roi_beta(x,y,z);
                            trial_roi_mean(1,trial_idx(b)) = nanmean(A);
                        end
                    end
                    roi_mean(i_sbj,:)=trial_roi_mean;

                end

            end
            current_con{r,1} = trials(t);
            current_con{r,2} = roi_masks{r};
            current_con{r,3} = roi_mean;
            T.(roi_masks{r}) = roi_mean;
        end
        conditions_roi = vertcat(conditions_roi, current_con);
        writetable(T, fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.xlsx')), 'Sheet', strcat('Trial_', num2str(t)));
    end
    save(fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.mat')),"current_con");


    save(fullfile(path_condition,strcat(phase,'_conditions_roi.mat')),"conditions_roi");

end




%% trial
% 
% if strcmp(cond_or_trial,'trial')
%     con_name=trial;
%     beta=trial_beta;
%     trials=1:32;
% 
%     n_tri=arrayfun(@(x) sprintf('Trial_%d',x),trials,'UniformOutput',false);
% 
%     conditions_roi={};
% 
%     for con_i=1:numel(con_name)
%         disp([cond_or_trial ' : ' con_name{con_i}])
%         current_con={}; % 여기서 초기화하기!!!
%         path_condition = fullfile(path_root,'ROI_result',(main_or_ODT),(phase),(mov_ver),(cond_or_trial));
%         mkdir(path_condition);
% 
%         for t=trials
%             current_trial={}; % 여기서 초기화하기!!!
%             % T = array2table(nan(numel(sbj_id_list), numel(trials)),'VariableNames', n_tri);
%             % T = array2table(nan(numel(sbj_id_list), numel(trials)));
%             roi_mean=[];
%             trial_roi_mean=zeros(1,32);
% 
%             for r=1:numel(which_roi)
%                 current_roi = fullfile(roi_path,roi_files(which_roi(r)).name);
%                 roi_coord = mask2coord(current_roi);
%                 for i_sbj = 1: numel(sbj_id_list)
%                     sbj_beta = beta.(con_name{con_i}){i_sbj,1};
%                     trial_idx = beta.(con_name{con_i}){i_sbj,2};
%                     if strcmp(con_name{con_i},'total')
%                         trial_idx=trials; end
% 
%                     for b=1:numel(sbj_beta)
%                         roi_beta = niftiread(sbj_beta{b});
% 
%                         A=[];
%                         for i=1:size(roi_coord,1)
%                             x=roi_coord(i,1); y=roi_coord(i,2); z=roi_coord(i,3);
%                             A(i)=roi_beta(x,y,z);
%                             trial_roi_mean(1,trial_idx(b)) = nanmean(A);
%                         end
% 
%                         roi_mean(i_sbj,:)=trial_roi_mean;
% 
%                     end
% 
% 
%                     % current_con{r,1} = roi_masks{r};
%                     % current_con{r,2} = {roi_mean};
%                     T= array2table(roi_mean,'VariableNames', n_tri);
%                 end
%                 % conditions_roi = vertcat(conditions_roi, current_con);
%                 writetable(T, fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.xlsx')), 'Sheet', roi_masks{r});
%             end
%         end
%     end
%     save(fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.mat')),"current_con");
% 
% 
%     save(fullfile(path_condition,strcat(phase,'_conditions_roi.mat')),"conditions_roi");
% 
% end
% 




%% new
if strcmp(cond_or_trial,'trial')
    con_name=trial;
    beta=trial_beta;
    trials=1:32;

    n_tri=arrayfun(@(x) sprintf('Trial_%d',x),trials,'UniformOutput',false);

    conditions_roi={};

    for con_i=1:numel(con_name)
        disp([cond_or_trial ' : ' con_name{con_i}])
        current_con={}; % 여기서 초기화하기!!!
        path_condition = fullfile(path_root,'ROI_result',(main_or_ODT),(phase),(mov_ver),(cond_or_trial));
        mkdir(path_condition);
        for r=1:numel(which_roi)
            current_roi = fullfile(roi_path,roi_files(which_roi(r)).name);
            roi_coord = mask2coord(current_roi);

            for t=trials
                current_trial={}; % 여기서 초기화하기!!!
                % T = array2table(nan(numel(sbj_id_list), numel(trials)),'VariableNames', n_tri);
                % T = array2table(nan(numel(sbj_id_list), numel(trials)));
                roi_mean=[];
                trial_roi_mean=zeros(1,32);

                for i_sbj = 1: numel(sbj_id_list)
                    sbj_beta = beta.(con_name{con_i}){i_sbj,1};
                    trial_idx = beta.(con_name{con_i}){i_sbj,2};
                    if strcmp(con_name{con_i},'total')
                        trial_idx=trials; end

                    for b=1:numel(sbj_beta)
                        roi_beta = niftiread(sbj_beta{b});

                        A=[];
                        for i=1:size(roi_coord,1)
                            x=roi_coord(i,1); y=roi_coord(i,2); z=roi_coord(i,3);
                            A(i)=roi_beta(x,y,z);
                            trial_roi_mean(1,trial_idx(b)) = nanmean(A);
                        end

                        roi_mean(i_sbj,:)=trial_roi_mean;

                    end

                    % current_con{r,1} = roi_masks{r};
                    % current_con{r,2} = {roi_mean};
                    T= array2table(roi_mean,'VariableNames', n_tri);
                end
                % conditions_roi = vertcat(conditions_roi, current_con);
                writetable(T, fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.xlsx')), 'Sheet', roi_masks{r});
            end
        end
        save(fullfile(path_condition,strcat(phase,'_',con_name{con_i},'.mat')),"current_con");
        save(fullfile(path_condition,strcat(phase,'_conditions_roi.mat')),"conditions_roi");

    end
end
