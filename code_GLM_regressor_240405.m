clear; clc; close;
path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\event_tsv_files';
% cd(path);

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

    curr_T = readtable(fullfile(path, 'event_table_MR.xlsx'),'Sheet',c_sbj);

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


    reg_for_glm{1,i} = struct('n_sbj',n_sbj,'main',main,'ODT_detail',ODT,'trial_detail',...
        struct('trial_hit',trial_hit,'trial_corr_rej',trial_corr_rej,'trial_miss',trial_miss,'trial_false',trial_false,'trial_to',trial_to));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make better regressor
% set path
% root_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';cd(root_path)
% load('reg_for_glm_mov')
%%
% merge regressors


for i = 1:length(sbj_id_list)
    %% main
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

    %% ODT
    O_onsets = reshape(reg_for_glm{1,i}.ODT_detail.regress_onset',[],1);
    O_durations = reshape(reg_for_glm{1,i}.ODT_detail.regress_duration',[],1);
    ODT_names = {'pre_forest','pre_city','post_forest','post_city','pre_target','post_target'}';
    ODT_onsets = {sort([O_onsets{1} O_onsets{2}]); sort([O_onsets{3} O_onsets{4}]); sort([O_onsets{5} O_onsets{6}]); sort([O_onsets{7} O_onsets{8}]); ...
        O_onsets{9}; O_onsets{10}};
    ODT_durations = {4;4;4;4;4;4};

    %% make MAIN suffix
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

    %% make ODT suffix
    pre_forest=repmat("pre_F_1",1,numel(ODT_onsets{1}));
    post_forest=repmat("post_F_1",1,numel(ODT_onsets{3}));
    pre_city=repmat("pre_C_1",1,numel(v2_onsets{2}));
    post_city=repmat("post_C_1",1,numel(v2_onsets{4}));

    pre_forest(~ismember(ODT_onsets{1},O_onsets{1})) = "pre_F_2";
    post_forest(~ismember(ODT_onsets{3},O_onsets{5})) = "post_F_2";
    pre_city(~ismember(ODT_onsets{2},O_onsets{3})) = "pre_C_2";
    post_city(~ismember(ODT_onsets{4},O_onsets{7})) = "post_C_2";

    pre_target=repmat("pre_target",1,numel(ODT_onsets{5}));
    post_target=repmat("post_target",1,numel(ODT_onsets{6}));

    ODT_suffix = {pre_forest;pre_city;post_forest;post_city;pre_target;post_target};

%% 0405 6am %%%%%%%%%%%%%%%


keys = {'o_hit', 'o_miss', 'o_false', 'o_corr_rej', 'o_timeout','c_hit', 'c_miss', ...
    'c_false', 'c_corr_rej', 'c_timeout', 'f_f', 'c_f', 'f_r', 'c_r'};
onset_val = {onsets{1}, onsets{4}, onsets{2}, onsets{3}, onsets{5}, onsets{6}, onsets{9}, ...
           onsets{7}, onsets{8}, onsets{10}, onsets{16}, onsets{17}, onsets{18}, onsets{19}};
suff_val = {repelem({'hit'}, numel(onsets{1})), repelem({'miss'}, numel(onsets{2})), repelem({'false'}, numel(onsets{3})), ...
          repelem({'corr_rej'}, numel(onsets{4})), repelem({'timeout'}, numel(onsets{5})), repelem({'hit'}, numel(onsets{6})), ...
          repelem({'miss'}, numel(onsets{7})), repelem({'false'}, numel(onsets{8})), repelem({'corr_rej'}, numel(onsets{9})), ...
          repelem({'timeout'}, numel(onsets{10})), repelem({'f_f'}, numel(onsets{16})), repelem({'c_f'}, numel(onsets{17})), ...
          repelem({'f_r'}, numel(onsets{18})), repelem({'c_r'}, numel(onsets{19}))};
mo = containers.Map(keys, onset_val);
ms = containers.Map(keys, suff_val);


mmo_key = {'corr';'incorr';'ccorr';'cincorr';'f_f';'c_f';'f_r';'c_r';'forest';'city'};
mmo_onset = {[mo('o_hit') mo('o_corr_rej')]; [mo('o_miss') mo('o_false') mo('o_timeout')]; [mo('c_hit') mo('c_corr_rej')]; [mo('c_miss') mo('c_false') mo('c_timeout')]; ...
        mo('f_f'); mo('c_f'); mo('f_r'); mo('c_r');[mo('f_f') mo('f_r')]; [mo('c_f') mo('c_r')]};
mms_sufff = {[ms('o_hit') ms('o_corr_rej')]; [ms('o_miss') ms('o_false') ms('o_timeout')]; [ms('c_hit') ms('c_corr_rej')]; [ms('c_miss') ms('c_false') ms('c_timeout')]; ...
        ms('f_f'); ms('c_f'); ms('f_r'); ms('c_r');[ms('f_f') ms('f_r')]; [ms('c_f') ms('c_r')]};
mmo = containers.Map(mmo_key, mmo_onset);
mms=containers.Map(mmo_key, mms_sufff);

total_key = {'obj_show'; 'choice';'mov1'; 'mov2'};
total_ons = {[mmo('corr') mmo('incorr')]; [mmo('ccorr') mmo('cincorr')]; [mmo('f_f') mmo('c_f') mmo('f_r') mmo('c_r')]; [mmo('forest') mmo('city')]};
total_sufff = {[mms('corr') mms('incorr')]; [mms('ccorr') mms('cincorr')]; [mms('f_f') mms('c_f') mms('f_r') mms('c_r')]; [mms('forest') mms('city')]};

total = containers.Map(total_key,total_ons);
ts=containers.Map(total_key,total_sufff);

% v1_suffix = {'corr';'incorr';'ccorr';'cincorr';'f_f';'c_f';'f_r';'c_r'};
% v2_suffix = {'corr';'incorr';'ccorr';'cincorr';'forest';'city'};


[sort_obj , sort_obj_idx] = sort(total('obj_show'));

total_obj_show = ts('obj_show');
obj_suffix=cellfun(@(x,y) sprintf('%s_%.2d',x,y),total_obj_show(sort_obj_idx), num2cell(1:32), 'UniformOutput', false);
[~, revert_idx] = sort(sort_obj_idx);
obj_suffix = obj_suffix(revert_idx);


[sort_choi , sort_choi_idx] = sort(total('choice'));

total_choice = ts('choice');
choic_suffix=cellfun(@(x,y) sprintf('%s_%.2d',x,y),total_choice(sort_choi_idx), num2cell(1:32), 'UniformOutput', false);
[~, revert_ch_idx] = sort(sort_choi_idx);
choic_suffix = choic_suffix(revert_ch_idx);



corr=obj_suffix(1:numel(mms('corr')));
incorr=obj_suffix(numel(mms('corr'))+1 : end);
ccorr=obj_suffix(1:numel(mms('ccorr')));
cincorr=obj_suffix(numel(mms('ccorr'))+1 : end);


v1_suffix = {corr;incorr;ccorr;cincorr;ms('f_f'); ms('c_f'); ms('f_r'); ms('c_r')};
v2_suffix = {corr;incorr;ccorr;cincorr;[ms('f_f') ms('f_r')]; [ms('c_f') ms('c_r')]};



    %% create reg_for_glm
    reg_for_glm{1, i}.ODT.regress_name = ODT_names;
    reg_for_glm{1, i}.ODT.regress_onset = ODT_onsets;
    reg_for_glm{1, i}.ODT.regress_duration = ODT_durations;
    reg_for_glm{1, i}.ODT.regress_suffix = ODT_suffix;
    reg_for_glm_ver1{1, i} = reg_for_glm{1, i}; reg_for_glm_ver2{1, i} = reg_for_glm{1, i};

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
save(fullfile('regressors_GLM_0513'),"reg_for_glm_ver1","reg_for_glm_ver2","sbj_id_list");
