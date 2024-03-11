clear; clc; close;
% addpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\new\data_bhv_log_table\GLM')
% load('event_regressor_all2');
% regressor= behav_regressor_glm{1,4}{1,1};
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
%     table_TR= readtable('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_log_table\total\event_TR.xlsx');
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
%
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
main_regress_name=cell(4,5);main_regress_duration=cell(4,5);
% name
main_regress_name(1,:) = {'obj_show_hit','obj_show_miss','obj_show_false','obj_show_corr_rej','obj_show_timeout'};
main_regress_name(2,:) = {'choice_hit','choice_miss','choice_false','choice_corr_rej', 'choice_timeout'};
main_regress_name(3,:) = {'obj_ITI_hit','obj_ITI_miss','obj_ITI_false','obj_ITI_corr_rej','obj_ITI_timeout'};
main_regress_name(4,:) = {'moving_forest', 'moving_city','tunnel_forest','tunnel_city',NaN};
% duration
main_regress_duration(1,:) = {4,4,4,4,4};
main_regress_duration(2,:) = {2,2,2,2,2};
main_regress_duration(3,:) = {2,2,2,2,2};
main_regress_duration(4,:) = {4,4,2,2,NaN};

%% OD task: regressor names cell
ODT_regress_name=cell(3,4);ODT_regress_duration=cell(3,4);
%name
ODT_regress_name(1,:) = {'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2'};
ODT_regress_name(2,:) = {'post_forest_1','post_forest_2','post_city_1','post_city_2'};
ODT_regress_name(3,:) = {'pre_all','post_all','fixation','target'};
%duration
ODT_regress_duration(1,:) = {4,4,4,4};
ODT_regress_duration(2,:) = {4,4,4,4};
ODT_regress_duration(3,:) = {120,120,NaN,4};

%%
reg_for_glm={};ODT_regress_onset=cell(3,4);
for i = 1: numel(sbj_id_list)
    n_sbj = sbj_id_list(i);
    c_sbj = sprintf('sub-%.2d',n_sbj);

    curr_T = readtable('event_table_MR.xlsx','Sheet',c_sbj);
    % sub-04의 post ODT에서 5번이 2번나오고 6번이 4번 나온 것으로 찍혀있음(24.03.22 5pm기준)
    % 우선 맨 마지막 trial을 6번이 아닌 5번으로 바꿔줘서 개수 맞춰줌. 저녁에 재민선배가 확인해주신다고 함
    if i==2
        curr_T.Obj_ID(curr_T.Trial==-30)=5;
    end

    % scan
    main_scan_num=[floor(curr_T.MR_start(curr_T.Trial==1)/TR), floor(curr_T.MR_end(curr_T.Trial==32)/TR)];
    pre_scan_num=[floor(curr_T.MR_start(curr_T.Trial==-1)/TR), floor(curr_T.MR_end(curr_T.Trial==-15)/TR)];
    post_scan_num=[floor(curr_T.MR_start(curr_T.Trial==-16)/TR), floor(curr_T.MR_end(curr_T.Trial==-30)/TR)];

    %% onset
    % associated object
    obj_f = unique(curr_T.Obj_ID(curr_T.Context_Num ==1 & curr_T.Association == 1));
    obj_c = unique(curr_T.Obj_ID(curr_T.Context_Num ==2 & curr_T.Association == 1));

    f1 = curr_T.ObjOn(curr_T.Obj_ID == min(obj_f));
    f2 = curr_T.ObjOn(curr_T.Obj_ID == max(obj_f));
    c1 = curr_T.ObjOn(curr_T.Obj_ID == min(obj_c));
    c2 = curr_T.ObjOn(curr_T.Obj_ID == max(obj_c));

    % ODT
    for odt_row=1:2
        if odt_row == 1; obj_idx=1:7;
        elseif odt_row == 2; obj_idx = 8:14;
        end
        ODT_regress_onset{odt_row,1} =f1(obj_idx);
        ODT_regress_onset{odt_row,2} =f2(obj_idx);
        ODT_regress_onset{odt_row,3} =c1(obj_idx);
        ODT_regress_onset{odt_row,4} =c2(obj_idx);
    end
    ODT_regress_onset{3,1} = curr_T.ObjOn(curr_T.Lap==0);
    ODT_regress_onset{3,2} = curr_T.ObjOn(curr_T.Lap==9);
    ODT_regress_onset{3,3} = curr_T.ObjOff(curr_T.Trial < 0);
    ODT_regress_onset{3,4} = curr_T.ObjOn(curr_T.Obj_ID==12);


    % onset
    onsets=[];main_regress_onset=cell(4,5);fix=[];
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
                    % Correct_rejection
                elseif curr_T.Association(row) ==0
                    main_regress_onset{reg_row,4}(end+1) = onsets(row);
                end
            elseif curr_T.Correct_Num(row) ==0
                % miss
                if curr_T.Association(row) ==1
                    main_regress_onset{reg_row,2}(end+1) = onsets(row);
                    % false alarm
                elseif curr_T.Association(row) ==0
                    main_regress_onset{reg_row,3}(end+1) = onsets(row);
                end
                % time out
            elseif curr_T.Correct_Num(row) ==2
                main_regress_onset{reg_row,5}(end+1) = onsets(row);
            end
        end
        if curr_T.Trial(row) > 0
            if curr_T.Context_Num(row) ==1 % forest
                col1=1; col2=3;
            elseif curr_T.Context_Num(row) ==2 % city
                col1=2; col2=4;
            end
            % Moving
            main_regress_onset{4,col1}(end+1) = curr_T.NavStart(row);
            % tunnel
            if curr_T.Lap_Trial(row) ==4
                main_regress_onset{4,col1}(end+1) = curr_T.ITIEnd(row);
                main_regress_onset{4,col2}(end+1) = curr_T.NavEnd(row);
            end
            %% ODT
            if curr_T.Trial(row) < 0 && row ~= 1
                % ODT duration
                fix = [fix; curr_T.ObjOn(row) - curr_T.ObjOff(row-1)];
            end

        end
        ODT_regress_duration{3,3} = fix;


        main=struct('scan_num',main_scan_num,'regress_name',{main_regress_name},...
            'regress_onset',{main_regress_onset},'regress_duration',{main_regress_duration});
        ODT= struct('pre_scan_num',pre_scan_num,'post_scan_num',post_scan_num,'regress_name',{ODT_regress_name},...
            'regress_onset',{ODT_regress_onset},'regress_duration',{ODT_regress_duration});


        reg_for_glm{1,i} = struct('n_sbj',n_sbj,'main',main,'ODT',ODT);
    end

end
save("reg_for_glm","reg_for_glm")
save("Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code\reg_for_glm","reg_for_glm")








    %
    % main_regress_onset{1,:} =
    % main_regress_onset{2,:} =
    % main_regress_onset{3,:} =
    % main_regress_onset{4,:} =


    % main=struct('scan_num',main_scan_num,'regress_name',{main_regress_name},...
    %     'regress_onset',{main_regress_onset},'regress_duration',{main_regress_duration});
    % preODT= struct('scan_num',pre_scan_num,'regress_name',{pre_regress_name},...
    %     'regress_onset',{pre_regress_onset},'regress_duration',{pre_regress_duration});
    %
    % postODT= struct('scan_num',post_scan_num,'regress_name',{post_regress_name},...
    %     'regress_onset',{post_regress_onset},'regress_duration',{post_regress_duration});
    %
    %
    %
    % reg_for_glm{1,n_sbj} = struct('n_sbj',n_sbj,'main',main,'preODT',preODT,'postODT',postODT);
    %
    %
    %
    %









