function [all_sbj_events_temp,num_sbj_events_temp,fig,sbj_info_file_temp]=func_bhv_logparsing(path,task_name,sbj_i,is_save_output,is_open_plot,sbj_info_file_temp)
all_sbj_events_temp=[];num_sbj_events_temp = [];

% set defaults
if ~exist('path','var'),  error('path is missing!'); end
if ~exist('task_name','var'), task_name='OCAT'; disp('Use default: task_name : OCAT'); end
if ~exist('sbj_i','var'),  error('sbj_i is missing!'); end
if ~exist('is_save_output','var'), error('is_save_output is missing!'); end
if ~exist('is_open_plot','var'), error('is_open_plot is missing!'); end


%% import the log file
flag_log = dir(fullfile(path{1},'*Behavior.csv'));
file_list = arrayfun(@(x) readtable(fullfile(x.folder, x.name)), flag_log, 'uni',0);

%% For Loop!!
    c_sbj = strcat('sub-', num2str(sbj_i, '%02.f'));

%% make output directory
   path_out = {};
   path_out{end+1} = fullfile(path{2},'individual',c_sbj);
   path_out{end+1} = fullfile(path{2},'total');
   path_out{end+1} = fullfile(path{2},'GLM');
   path_out{end+1} = fullfile(path{3},'performance');
   path_out{end+1} = fullfile(path{4},c_sbj,'func');
   path_out{end+1} = fullfile(path{5},'Responses');
     
   
      for i=1:length(path_out)
         mkdir(path_out{i});
      end
  

    %% remove time from sbj events table
    sbj_events = file_list{sbj_i};
    sbj_events(contains(sbj_events.Var1(:),'time'),:) = [];


    %% save TR log (based on the first MR imaging time)
    event_MR = find(contains(sbj_events.Var1(:),'MR'));
    first_event_MR = sbj_events.Var2(event_MR(1));

    idx = sbj_events.Var2 ~= fix(sbj_events.Var2);
    sbj_events.Var2(idx) = sbj_events.Var2(idx) - first_event_MR;

    event_TR = sbj_events(event_MR,:);
    if is_save_output == 1; writetable(event_TR,[path_out{1} '\' c_sbj '_event_TR.csv']);end
    % remove MR row
    sbj_events(event_MR,:) = [];


    %% save PV task log
    PV_boundary = find(contains(sbj_events.Var1(:),'OCP'));
    event_pre_PV = sbj_events(PV_boundary(1):PV_boundary(2),:);
    event_post_PV = sbj_events(PV_boundary(3):PV_boundary(4),:);
    if is_save_output == 1
    writetable(event_pre_PV,[path_out{1} '\' c_sbj '_event_pre_PV.csv']);
    writetable(event_post_PV,[path_out{1} '\' c_sbj '_event_post_PV.csv']);
    end
    
     %% DMTS task 중 target에 대한 RT
    % pre-PV와 post-PV 이벤트를 합칩니다.
    sbj_DMTS = [event_pre_PV; event_post_PV];

    % Var4 값이 '12'이고 그 다음 행의 Var1이 'ButtonA'인 경우를 찾습니다.
    DMTS_RT=[];
    for idx=1:height(sbj_DMTS)-1
        if sbj_DMTS.Var4(idx)==12 && strcmp(sbj_DMTS.Var1(idx+1), 'ButtonA')
            DMTS_RT = [DMTS_RT; sbj_DMTS.Var2(idx + 1) - sbj_DMTS.Var2(idx)];
        end
    end
    RT_mean = mean(DMTS_RT);

    % 'DMTS_RT'라는 새로운 변수를 sbj_info_file 테이블에 추가하고, 그 변수의 sbj_i행에 평균값을 넣습니다.
    sbj_info_file_temp.DMTS_RT(sbj_i) = RT_mean;

    %% DMTS event onset --> event table에 추가
    % pre는 Lap== 0; trial -1 ~ -15 // post는 Lap==9, trial -16 ~ -30
    pre_DMTS=event_pre_PV(strcmp(event_pre_PV.Var1(:),'PreObjOn'),:);
    pre_DMTS.Var1=zeros(height(pre_DMTS),1); % Lap ==0
    pre_DMTS.Var3=[-1:-1:-1*height(pre_DMTS)]'; % trial== -1:-15

        
    post_DMTS=event_post_PV(strcmp(event_post_PV.Var1(:),'PreObjOn'),:);
    post_DMTS.Var1=9*ones(height(post_DMTS),1); % Lap ==9
    post_DMTS.Var3=[-16:-1:(-1*(height(post_DMTS)+15))]'; % trial== -1:-15
        
    %% object #4,5,6,7 only!
    pre_DMTS_obj = pre_DMTS((pre_DMTS.Var4(:)~=12),:);
    post_DMTS_obj = post_DMTS((post_DMTS.Var4(:)~=12),:);
    
    
%% Define variable names
var_name = ["Lap", "TrialStart", "Trial","Lap_Trial", "Context_txt", "Context_Num", ...
"Direction", "Location","Association", "Obj_ID", "ObjOn","ChoiceOn", "Choice_Num", ...
"Choice_txt","Correct_Num","Correct_txt","isTimeout" "RT", "ObjOff", "TrialEnd"];

event_struct = struct;
for v = 1:length(var_name)
    if contains(var_name{v}, '_txt')
        event_struct.(var_name{v}) = {};
    else
        event_struct.(var_name{v}) = [];
    end
end


for i = 6:10
    type_log = find(contains(sbj_events.Var3(:),'Type'));
    type_log = num2str(sbj_events.Var4(type_log), i-5) - '0';
end
    for j = 6:10
    event_struct.(var_name{j}) = type_log(:,j-5);
    end
event_struct.Lap_Trial(1,1)=1;



%% parsing events
event_name = string(sbj_events.Var1(:));
% Lap & Trial
time_lap_start = sbj_events(event_name=="LapStart",[2 4]);
for i=1:height(sbj_events)
    if event_name(i)=="TrialStart"
        event_struct.TrialStart(end+1,1) = sbj_events{i,2};
    elseif event_name(i)=="Trial"
        event_struct.Trial(end+1,1) = sbj_events{i,2};
        event_struct.Lap_Trial(end+1,1) = mod(length(event_struct.Trial), 4)+1;
        lap_idx = find(event_struct.TrialStart(end) > time_lap_start.Var2(:), 1, 'last');
        if ~isempty(lap_idx) && event_struct.TrialStart(end) > time_lap_start.Var2(lap_idx)
         event_struct.Lap(end+1,1) = time_lap_start.Var4(lap_idx);
        end
       

   % Choice
    else
        if contains(event_name(i), "Choice"+("A"|"B"))
        event_struct.Choice_txt{end+1,1} = string(extractAfter(event_name(i),"Choice"));
   % correctness, RT, isTimeout
        elseif event_name(i)=="Decision"
        %correct Timeout error
                if sbj_events{i,4} > 1.5 && sbj_events{i,2} ~= 2
                    sbj_events{i,2} = 2;
                end

             %correctness
             event_struct.Correct_Num(end+1,1) = sbj_events{i,2};
             %RT
             event_struct.RT(end+1,1) = sbj_events{i,4};

            %txt, timeout
           switch sbj_events{i,2}
               case 0
                event_struct.Correct_txt = [event_struct.Correct_txt; "Incorrect"];
                event_struct.isTimeout = [event_struct.isTimeout; 0];
               case 1
                event_struct.Correct_txt = [event_struct.Correct_txt; "Correct"];
                event_struct.isTimeout = [event_struct.isTimeout; 0];
               case 2
                event_struct.Correct_txt = [event_struct.Correct_txt; "TimeOut"];
                event_struct.isTimeout = [event_struct.isTimeout; 1];
           end
    
        %rest of fields    
        elseif ismember(event_name(i),var_name) && ~strcmp(event_name(i), 'Trial')
             event_struct.(event_name(i)){end+1,1} = sbj_events{i,2};
        end
    end

    %choice missing (choiceOn과 objOff 이벤트 개수와 choice 개수가 다를 경우 missing 삽입)
    if length(event_struct.ChoiceOn) == length(event_struct.ObjOff) && length(event_struct.ChoiceOn) ~= length(event_struct.Choice_txt)
    event_struct.Choice_txt{end+1,1} = missing;
    end
end

% Choice txt to Number
event_struct.Choice_Num = NaN(height(event_struct.Choice_txt), 1);  % Initialize Choice_Num with NaN
event_struct.Choice_Num = cellfun(@(x) find(strcmpi(x, {'A', 'B'}), 1), event_struct.Choice_txt, 'UniformOutput', false);

% Context Num to txt
event_struct.Context_txt = replace(string(event_struct.Context_Num), ["1", "2"], ["F", "C"]);

%% Making a Table
event_struct = orderfields(event_struct,var_name);
event_struct.Lap_Trial(end)=[];
event_table = struct2table(event_struct);

%% put pre-DMTS into event table

pre_temp_table = array2table(nan(height(pre_DMTS),length(var_name)));pre_temp_table.Properties.VariableNames = cellstr(var_name);
pre_temp_table(:,[1,2,3,10,11]) = pre_DMTS(:,[1,2,3,4,2]);

post_temp_table = array2table(nan(height(pre_DMTS),length(var_name)));post_temp_table.Properties.VariableNames = cellstr(var_name);
post_temp_table(:,[1,2,3,10,11]) = post_DMTS(:,[1,2,3,4,2]);

% ObjID가 어느 context에 속하는지 여부 찾아서 Context_Num, Context_txt, Association 열 채우기

tables = {pre_temp_table, post_temp_table};

for t = 1:length(tables)
    current_table = tables{t};
    current_table.Context_txt=string(current_table.Context_txt);
    unique_values = unique(current_table.Obj_ID);

    % 각 고유값에 대해 첫 번째 일치하는 행 찾기
    for value = unique_values'
        % event_table.Var10이 value와 일치하고 event_table.Var9가 1인 첫 번째 행 찾기
        matching_row = find(event_table.Obj_ID == value & event_table.Association == 1, 1, 'first');
        if ~isempty(matching_row)
            current_table.Context_txt(current_table.Obj_ID == value) = string(event_table.Context_txt(matching_row));
            current_table.Context_Num(current_table.Obj_ID == value) = event_table.Context_Num(matching_row);
            current_table.Association(current_table.Obj_ID == value) = event_table.Association(matching_row);
        end
    end
    tables{t} = current_table;
end

% 결과를 pre_DMTS_table과 post_DMTS_table에 저장
pre_temp_table = tables{1};
post_temp_table = tables{2};


var_name_num = ["ObjOn", "ChoiceOn", "Choice_Num", "ObjOff","TrialEnd"];

for i = 1:length(var_name_num)
   event_table.(var_name_num(i)) = cell2mat(event_table.(var_name_num(i)));
end
event_table.Choice_txt = string(event_table.Choice_txt);

event_table = vertcat(pre_temp_table, event_table, post_temp_table);

%% save the table
if is_save_output == 1
%individual
writetable(event_table,[path_out{1} '\' c_sbj '_event_table.csv']);
% total table
    %TR
    writetable(event_TR,[path_out{2} '\event_TR.xlsx'],'Sheet',c_sbj);
    %PV
    writetable(event_pre_PV,[path_out{2} '\event_pre_PV.xlsx'],'Sheet',c_sbj);
    writetable(event_post_PV,[path_out{2} '\event_post_PV.xlsx'],'Sheet',c_sbj);
    %events
    writetable(event_table,[path_out{2} '\event_table.xlsx'],'Sheet',c_sbj);
    
% bids .tsv
writetable(event_table,[path_out{5} '\' c_sbj '_task-ocat_run-01_events.tsv']);
end

% correct_regressor; object가 켜진 시간, for GLM
corr_match = cell2mat(event_struct.ObjOn(event_struct.Correct_Num == 1 & event_struct.Association ==1));
corr_nonmatch = cell2mat(event_struct.ObjOn(event_struct.Correct_Num == 1 & event_struct.Association ==0));

incorr_match = cell2mat(event_struct.ObjOn(event_struct.Correct_Num ~= 1 & event_struct.Association ==1));
incorr_nonmatch = cell2mat(event_struct.ObjOn(event_struct.Correct_Num ~= 1 & event_struct.Association ==0)); 

%% Context-Object 조합
    combi_C = event_struct.Context_txt(event_struct.Association ==1);
    combi_O = event_struct.Obj_ID(event_struct.Association ==1);

    combi_FFCC = [sort(combi_O(find(combi_C=="F",2,"first"))); sort(combi_O(find(combi_C=="C",2,"first")))]';
    ffcc_array = num2cell(combi_FFCC);
    combi_FFCC = strjoin(string(combi_FFCC));
    disp(['combi_FFCC: ', combi_FFCC]);
    sbj_info_file_temp.Combi_FFCC{sbj_i} = combi_FFCC;
    
    onset_pre_DMTS_forest = table2array(pre_DMTS_obj((pre_DMTS_obj.Var4(:)==ffcc_array{1} | pre_DMTS_obj.Var4(:)==ffcc_array{2}),2));
    onset_pre_DMTS_city = table2array(pre_DMTS_obj((pre_DMTS_obj.Var4(:)==ffcc_array{3} | pre_DMTS_obj.Var4(:)==ffcc_array{4}),2));
    onset_post_DMTS_forest = table2array(post_DMTS_obj((post_DMTS_obj.Var4(:)==ffcc_array{1} | post_DMTS_obj.Var4(:)==ffcc_array{2}),2));
    onset_post_DMTS_city = table2array(post_DMTS_obj((post_DMTS_obj.Var4(:)==ffcc_array{3} | post_DMTS_obj.Var4(:)==ffcc_array{4}),2));
    
    %% 모든 condition들을 cell array로 묶기
    names = {'corr_match', 'corr_nonmatch', 'incorr_match', 'incorr_nonmatch', 'onset_pre_DMTS_forest', 'onset_pre_DMTS_city', 'onset_post_DMTS_forest', 'onset_post_DMTS_city'};
    onsets = {corr_match, corr_nonmatch, incorr_match, incorr_nonmatch, onset_pre_DMTS_forest, onset_pre_DMTS_city, onset_post_DMTS_forest, onset_post_DMTS_city};
    durations = {4, 4, 4, 4, 4, 4, 4, 4}; % 모든 duration은 4초

    if is_save_output == 1
       save(fullfile(path_out{3},[c_sbj, '_multiple_conditions.mat']), 'names', 'onsets', 'durations');
    end
    
    %% learning curve 용 correct .mat file 제작 Responses -> 0 : incorrect trial + timeout, 1: correct trial
    Responses = event_struct.Correct_Num';
    Responses(Responses==2) = 0;     
   if is_save_output == 1
    save(fullfile(path_out{6},[c_sbj, '_Responses.mat']), 'Responses');
   end
   
    %% Movement Regressor 저장
%% Path containing data
% path for confounding factors
conf_file = dir(fullfile(path{4},'derivatives', c_sbj,'func', ['*task-', task_name, '*confounds*.tsv']));   % find the file

%% create "R" variable from movement_regressor matrix and save
    movereg_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};

    [data1, header1, ] = tsvread(fullfile(conf_file.folder, conf_file.name));
    
    trans_x=strmatch(movereg_names{1},header1,'exact');
    trans_y=strmatch(movereg_names{2},header1,'exact');
    trans_z=strmatch(movereg_names{3},header1,'exact');
    rot_x=strmatch(movereg_names{4},header1, 'exact');
    rot_y=strmatch(movereg_names{5},header1,'exact');
    rot_z=strmatch(movereg_names{6},header1,'exact');
    
    % remove the first row, 26-31 columns --> movement regressors    R_mov1 = fillmissing(R_mov1, ?nearest?);    R_mov1(~isfinite(R_mov1))=0;
    R = data1(2:end, [trans_x,trans_y,trans_z,rot_x,rot_y,rot_z]);
    R = fillmissing(R, 'nearest');
    R(~isfinite(R)) = 0;
    
    %R = data1(2:end, (end-5):end);  % remove the first row, 26-31 columns --> movement regressors
    if is_save_output == 1
    save(fullfile(path_out{3},[c_sbj, '_movereg.mat']), 'R');
    end    
    
    
%% Table for Analysis == event_table_NumOnly
session = repmat(sbj_i, 32,1);
event_numeric = struct;
event_numeric.session = session;

var_name_num = ["Lap", "Trial", "Context_Num", "Direction", "Location", "Association", "Obj_ID", "Choice_Num", "Correct_Num", "RT", "isTimeout"];
new_var_name = ["Lap", "Trial", "Context", "Direction", "Location", "Association", "Object", "Choice", "Correct", "RT", "isTimeout"];

for i = 1:length(var_name_num)
   event_numeric.(new_var_name(i)) = event_struct.(var_name_num(i));
end

event_numeric = struct2table(event_numeric);


% save /a subject
if is_save_output == 1
    writetable(event_numeric,[path_out{1} '\' c_sbj '_event_numeric_table.csv']);
end
% save /all subjects
num_sbj_events_temp = [num_sbj_events_temp;event_numeric];
all_sbj_events_temp = [all_sbj_events_temp;event_table];

%% %%%%%%% performance plot (accuracy, RT) %%%%%%%%%
    % Create a new figure for every 4 subjects   
        if is_open_plot == 1
            figure('Position', [1500 500 1000 600]);
        else
            figure('Position', [1500 500 1000 600], 'Visible', 'off');
        end

      %%% Set the current figure to f1
        hold on
        graph = func_perform_graph(event_table, c_sbj);
        f1=graph;
        % Save the plot for the current subject
        if is_save_output == 1
            saveas(f1, [path_out{4} '\individual\' c_sbj '_Performance.png']);
        end
        hold off
    fig{sbj_i}=f1;

end

    
    
    
    
