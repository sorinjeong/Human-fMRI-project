clear; clc; close;
path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\event_tsv_files';
cd(path);

% %%%%%%%% make events log table %%%%%%%%%%
%
% %% edit event table --> add Navigation period
% files=dir(fullfile(path, '*events.tsv'));
% for i=1:numel(files)
% T=readtable(files(i).name,"FileType","text",'Delimiter', '\t');
% T.TrialStart(T.Lap~=0 & T.Lap~=9)=T.TrialStart(T.Lap~=0 & T.Lap~=9)-4;
% T=renamevars(T,["TrialStart", "TrialEnd"],["NavStart", "ITIEnd"]);
% T.NavEnd=NaN(height(T),1);
% T.NavEnd(T.Lap_Trial==4)=T.ITIEnd(T.Lap_Trial==4)+4;
% writetable(T,'all_subjects_task-ocat_run-01_events.xlsx','Sheet',sprintf('sub-%.2d',i))
% end
%
%
%% align TR with events
% x="event_table_MR.xlsx";
%
% [~,sheets] = xlsfinfo(x);
%
% for s=1:numel(sheets)
%     table_events= readtable(x,'Sheet',sheets{s});
%     table_TR= readtable('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_log_table\total\event_TR.xlsx','Sheet',sheets{s});
%
%     table_events.MR_start=zeros(height(table_events),1);
%     table_events.MR_end=zeros(height(table_events),1);
%
%     for i= 1:height(table_events)
%         [~,idx1]=min(abs(table_TR.Var2 - table_events.NavStart(i)));
%         if table_events.Trial(i) <0
%             [~,idx2]=min(abs(table_TR.Var2 - table_events.ObjOff(i)));
%         else
%             [~,idx2]=min(abs(table_TR.Var2 - table_events.NavEnd(i)));
%         end
%         table_events.MR_start(i)=table_TR.Var2(idx1);
%         table_events.MR_end(i)=table_TR.Var2(idx2);
%
%     end
%
%     writetable(table_events, 'event_table_MR.xlsx','Sheet',sheets{s});
% end
% %
% %% add ODT object off onset
% for n_sub = 1:31
%     c_sub=sprintf('sub-%.2d',n_sub);
%
% A = readtable('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_log_table\total\event_pre_PV.xlsx', 'Sheet', c_sub);
% B = readtable('event_table_MR.xlsx', 'Sheet', c_sub);
% C = readtable('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_log_table\total\event_post_PV.xlsx', 'Sheet', c_sub);
%
%
% preoff = A.Var2(strcmp(A.Var1,'PreObjOff'));
% B.ObjOff(1:15) = preoff;
% postoff = C.Var2(strcmp(C.Var1,'PreObjOff'));
% B.ObjOff(48:62) = postoff;
%
% writetable(B, 'event_table_MR.xlsx', 'Sheet', c_sub);
% end

%% %%%%%%% make regressor cell array %%%%%%%%%%%%%%%%%%
% sbj_id_list = sbj_id_list(sbj_id_list~=7); % #7 제외
% save("sbj_id_list","sbj_id_list")


%% INPUT!!
TR = 2;

load('sbj_id_list');
%% main task: regressor names cell
main_regress_name={};main_regress_duration={};
% name
main_regress_name(1,:) = {'obj_show_hit','obj_show_miss','obj_show_false','obj_show_corr_rej','obj_show_timeout'};
main_regress_name(2,:) = {'choice_hit','choice_miss','choice_false','choice_corr_rej', 'choice_timeout'};
main_regress_name(3,:) = {'obj_ITI_hit','obj_ITI_miss','obj_ITI_false','obj_ITI_corr_rej','obj_ITI_timeout'};
main_regress_name(4,:) = {'moving_forest_first', 'moving_city_first','moving_forest_remainder','moving_city_remainder','tunnel'};
% duration
main_regress_duration(1,:) = {4,4,4,4,4};
main_regress_duration(2,:) = {2,2,2,2,2};
main_regress_duration(3,:) = {2,2,2,2,2};
main_regress_duration(4,:) = {4,4,4,4,2};

%% OD task: regressor names cell
ODT_regress_name=cell(3,4);ODT_regress_duration=cell(3,4);
%name
ODT_regress_name(1,:) = {'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2'};
ODT_regress_name(2,:) = {'post_forest_1','post_forest_2','post_city_1','post_city_2'};
ODT_regress_name(3,:) = {'pre_target','post_target','fixation',NaN};
%duration
% pre-ODT_all = 120s; post-ODT_all = 120s
ODT_regress_duration(1,:) = {4,4,4,4};
ODT_regress_duration(2,:) = {4,4,4,4};
ODT_regress_duration(3,:) = {4,4,4,NaN};

%%
reg_for_glm={};ODT_regress_onset=cell(3,4);
for i = 1: numel(sbj_id_list)
    n_sbj = sbj_id_list(i);
    c_sbj = sprintf('sub-%.2d',n_sbj);

    curr_T = readtable('event_table_MR.xlsx','Sheet',c_sbj);

    % onset after trim
    main_onset= curr_T.NavStart(curr_T.Trial==1);
    main_offset = curr_T.NavEnd(curr_T.Trial==32);
    pre_onset = curr_T.ObjOn(curr_T.Trial==-1);
    pre_offset = curr_T.ObjOff(curr_T.Trial==-15);
    post_onset = curr_T.ObjOn(curr_T.Trial==-16);
    post_offset = curr_T.ObjOff(curr_T.Trial==-30);

    % scan
    main_scan_num=[floor(main_onset/TR)+1, round(main_offset/TR)+1];
    pre_scan_num=[floor(pre_onset/TR)+1, round(pre_offset/TR)+1];
    post_scan_num=[floor(post_onset/TR)+1, round(post_offset/TR)+1];

    %% onset
    % associated object
    obj_f = unique(curr_T.Obj_ID(curr_T.Context_Num ==1 & curr_T.Association == 1));
    obj_c = unique(curr_T.Obj_ID(curr_T.Context_Num ==2 & curr_T.Association == 1));

    f1 = curr_T.ObjOn(curr_T.Obj_ID == min(obj_f) & curr_T.Trial < 0);
    f2 = curr_T.ObjOn(curr_T.Obj_ID == max(obj_f) & curr_T.Trial < 0);
    c1 = curr_T.ObjOn(curr_T.Obj_ID == min(obj_c) & curr_T.Trial < 0);
    c2 = curr_T.ObjOn(curr_T.Obj_ID == max(obj_c) & curr_T.Trial < 0);

    % ODT - pre 바로 다음에 post로 올 것 계산해서 onset 맞추기
    for odt_row=1:2
        if odt_row == 1; obj_idx=1:3; ons=pre_onset;
        elseif odt_row == 2; obj_idx = 4:6; ons=post_onset-(pre_offset+TR);
        end
        ODT_regress_onset{odt_row,1} =f1(obj_idx)-ons;
        ODT_regress_onset{odt_row,2} =f2(obj_idx)-ons;
        ODT_regress_onset{odt_row,3} =c1(obj_idx)-ons;
        ODT_regress_onset{odt_row,4} =c2(obj_idx)-ons;
    end
    target = curr_T.ObjOn(curr_T.Obj_ID==12);
    ODT_regress_onset{3,1} = target(1:3) - pre_onset;
    ODT_regress_onset{3,2} = target(4:6)  - post_onset+pre_offset+TR;

    ODT_regress_onset{3,3}= curr_T.ObjOff(curr_T.Trial < 0);
    ODT_regress_onset{3,3}(1:15) = ODT_regress_onset{3,3}(1:15) - pre_onset;
    ODT_regress_onset{3,3}(16:end) = ODT_regress_onset{3,3}(16:end) - post_onset+pre_offset+TR;

    ODT_regress_onset{3,4}=NaN;


    % onset
    onsets=[];main_regress_onset=cell(4,5);fix=[];
    trial_hit=[];trial_corr_rej=[];trial_miss=[];trial_false=[];trial_to=[];
  
   
    for row = 1:height(curr_T)
        %% main task
        for reg_row = 1:3
            if reg_row == 1; onsets = curr_T.ObjOn;
            elseif reg_row == 2; onsets = curr_T.ChoiceOn;
            elseif reg_row == 3; onsets = curr_T.ObjOff;
            end

            if curr_T.Correct_Num(row) ==1
                % Hit
                if curr_T.Association(row) ==1
                    main_regress_onset{reg_row,1}(end+1) = onsets(row);
                    if reg_row == 1;trial_hit(end+1) = curr_T.Trial(row);end
                    % Correct_rejection
                elseif curr_T.Association(row) ==0
                    main_regress_onset{reg_row,4}(end+1) = onsets(row);
                    if reg_row == 1;trial_corr_rej(end+1) = curr_T.Trial(row);end

                end
            elseif curr_T.Correct_Num(row) ==0
                % miss
                if curr_T.Association(row) ==1
                    main_regress_onset{reg_row,2}(end+1) = onsets(row);
                    if reg_row == 1;trial_miss(end+1) = curr_T.Trial(row);end

                    % false alarm
                elseif curr_T.Association(row) ==0
                    main_regress_onset{reg_row,3}(end+1) = onsets(row);
                    if reg_row == 1;trial_false(end+1) = curr_T.Trial(row);end

                end
                % time out
            elseif curr_T.Correct_Num(row) ==2
                main_regress_onset{reg_row,5}(end+1) = onsets(row);
                if reg_row == 1;trial_to(end+1) = curr_T.Trial(row);end

            end
        end
        if curr_T.Trial(row) > 0
            if curr_T.Context_Num(row) ==1 % forest
                col1=1; col2=3;
            elseif curr_T.Context_Num(row) ==2 % city
                col1=2; col2=4;
            end
            % Moving
            if curr_T.Lap_Trial(row) ==1
                main_regress_onset{4,col1}(end+1) = curr_T.NavStart(row);
            else
                main_regress_onset{4,col2}(end+1) = curr_T.NavStart(row);
            end
            % tunnel
            if curr_T.Lap_Trial(row) ==4
                main_regress_onset{4,col2}(end+1) = curr_T.ITIEnd(row);
                main_regress_onset{4,5}(end+1) = curr_T.NavEnd(row);
            end
        end
        %% ODT
        if curr_T.Trial(row) < 0 && row ~= 1
            % ODT duration
            fix = [fix; curr_T.ObjOn(row) - curr_T.ObjOff(row-1)];
        end

    end

    ODT_regress_duration{3,3} = fix;

    %% set onset to scan start time
    main_regress_onset = cellfun(@(x) x(~isnan(x))-main_onset, main_regress_onset, 'UniformOutput', false);

    %%

    main=struct('scan_num',main_scan_num,'regress_name',{main_regress_name},...
        'regress_onset',{main_regress_onset},'regress_duration',{main_regress_duration});
    ODT= struct('pre_scan_num',pre_scan_num,'post_scan_num',post_scan_num,'regress_name',{ODT_regress_name},...
        'regress_onset',{ODT_regress_onset},'regress_duration',{ODT_regress_duration});


    reg_for_glm{1,i} = struct('n_sbj',n_sbj,'main',main,'ODT',ODT,'trial_detail',...
        struct('trial_hit',trial_hit,'trial_corr_rej',trial_corr_rej,'trial_miss',trial_miss,'trial_false',trial_false,'trial_to',trial_to));
end


%% make better regressor
% set path
root_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';

%% 
% merge regressors
reg_for_glm_ver1 = reg_for_glm; reg_for_glm_ver2 = reg_for_glm;

for i = 1:length(sbj_id_list)

    % names = reshape(reg_for_glm{1,i}.((main_or_ODT)).regress_name',[],1);
    onsets = reshape(reg_for_glm{1,i}.main.regress_onset',[],1);
    durations = reshape(reg_for_glm{1,i}.main.regress_duration',[],1);

    v1_names = {'obj_show_corr','obj_show_incorr','choice_corr','choice_incorr','mov_forest_first','mov_city_first','mov_forest_remainder','mov_city_remainder'}';
    v1_onsets = {sort([onsets{1} onsets{4}]); sort([onsets{2} onsets{3} onsets{5}]); sort([onsets{6} onsets{9}]); sort([onsets{7} onsets{8} onsets{10}]); ...
        onsets{16}; onsets{17}; onsets{18}; onsets{19}};
    v1_durations = {4;4;2;2;4;4;4;4};

    v2_names = {'obj_show_corr','obj_show_incorr','choice_corr','choice_incorr','mov_forest','mov_city'}';
    v2_onsets = {sort([onsets{1} onsets{4}]); sort([onsets{2} onsets{3} onsets{5}]); sort([onsets{6} onsets{9}]); sort([onsets{7} onsets{8} onsets{10}]); ...
        sort([onsets{16} onsets{18}]); sort([onsets{17} onsets{19}])};
    v2_durations = {4;4;2;2;4;4};

%% make suffix
% corr=repmat("hit",1,numel(v1_onsets{1}));
% incorr=repmat("miss",1,numel(v1_onsets{2}));
forest=repmat("remainder",1,numel(v2_onsets{5}));
city=repmat("remainder",1,numel(v2_onsets{6}));

% for ii=[1 3]
%     for jj=[2 4]
%     corr(~ismember(v1_onsets{ii},onsets{1})) = "corr_rej";
%     incorr(ismember(v1_onsets{jj},onsets{3})) = "false";
%     incorr(ismember(v1_onsets{jj},onsets{5})) = "timeout";
%     end
% end
corr=cell(1,32);incorr=cell(1,32);
corr(reg_for_glm{1,i}.trial_detail.trial_hit) = {"hit"};
corr(reg_for_glm{1,i}.trial_detail.trial_corr_rej) = {"corr_rej"};
incorr(reg_for_glm{1,i}.trial_detail.trial_miss) = {"miss"};
incorr(reg_for_glm{1,i}.trial_detail.trial_false) = {"false"};
incorr(reg_for_glm{1,i}.trial_detail.trial_to) = {"timeout"};
corr(cellfun(@isempty,corr))=[];incorr(cellfun(@isempty,incorr))=[];

forest(~ismember(v2_onsets{5},onsets{18})) = "first";
city(~ismember(v2_onsets{6},onsets{19})) = "first";

f_f=repmat("f_f",1,numel(onsets{16}));
c_f=repmat("c_f",1,numel(onsets{17}));
f_r=repmat("f_r",1,numel(onsets{18}));
c_r=repmat("c_r",1,numel(onsets{19}));
v1_suffix = {corr;incorr;corr;incorr;f_f;c_f;f_r;c_r};
v2_suffix = {corr;incorr;corr;incorr;forest;city};

%% 



%% create reg_for_glm

reg_for_glm_ver1{1, i}.main.regress_name = v1_names;
reg_for_glm_ver1{1, i}.main.regress_onset = v1_onsets;
reg_for_glm_ver1{1, i}.main.regress_duration = v1_durations;
reg_for_glm_ver1{1, i}.main.regress_suffix = v1_suffix;

reg_for_glm_ver2{1, i}.main.regress_name = v2_names;
reg_for_glm_ver2{1, i}.main.regress_onset = v2_onsets;
reg_for_glm_ver2{1, i}.main.regress_duration = v2_durations;
reg_for_glm_ver2{1, i}.main.regress_suffix = v2_suffix;
end

%% save regressors
save(fullfile(root_path,'regressors_GLM_0402'),"reg_for_glm_ver1","reg_for_glm_ver2","sbj_id_list");
