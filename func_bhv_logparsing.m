function [all_sbj_events_temp,num_sbj_events_temp]=func_bhv_logparsing(path,task_name,sbj_i,is_save_output,is_open_plot)
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
sbj_info_file = readtable('../../OCAT_DIR/data/data_fmri_bids/participants.tsv','FileType','text');
sbj_info_file = removevars(sbj_info_file,["Weight","Size"]);

%% For Loop!!
    c_sbj = strcat('sub-', num2str(sbj_i, '%02.f'));

%% make output directory
   path_out = {};
   path_out{end+1} = fullfile(path{2},'individual',c_sbj);
   path_out{end+1} = fullfile(path{2},'total');
   path_out{end+1} = fullfile(path{2},'GLM');
   path_out{end+1} = fullfile(path{3},'performance');
     
    if ~exist(path{2},"dir")
      mkdir(path{2});mkdir(path{3});
      for i=2:length(path_out)
         mkdir(path_out{i});
      end
    end
    mkdir(path_out{1});

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
    sbj_info_file.DMTS_RT(sbj_i) = RT_mean;

    %% object #4,5,6,7 only!
    pre_DMTS_obj = event_pre_PV((event_pre_PV.Var4(:)==4 | event_pre_PV.Var4(:)==5 |event_pre_PV.Var4(:)==6 |event_pre_PV.Var4(:)==7),:);
    post_DMTS_obj = event_post_PV((event_post_PV.Var4(:)==4 | event_post_PV.Var4(:)==5 |event_post_PV.Var4(:)==6 |event_post_PV.Var4(:)==7),:);

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



%% parsing events
event_name = string(sbj_events.Var1(:));
% Lap & Trial
time_lap_start = sbj_events(event_name=="LapStart",[2 4]);
for i=1:height(sbj_events)
    if event_name(i)=="TrialStart"
        event_struct.TrialStart(end+1,1) = sbj_events{i,2};
    elseif event_name(i)=="Trial"
        event_struct.Trial(end+1,1) = sbj_events{i,2};
        event_struct.Lap_Trial(end+1,1) = mod(length(event_struct.Trial), 4) + 1;
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
event_table = struct2table(event_struct);

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
    sbj_info_file.Combi_FFCC{sbj_i} = combi_FFCC;
    
    onset_pre_DMTS_forest = table2array(pre_DMTS_obj((pre_DMTS_obj.Var4(:)==ffcc_array{1} | pre_DMTS_obj.Var4(:)==ffcc_array{2}),2));
    onset_pre_DMTS_city = table2array(pre_DMTS_obj((pre_DMTS_obj.Var4(:)==ffcc_array{3} | pre_DMTS_obj.Var4(:)==ffcc_array{4}),2));
    onset_post_DMTS_forest = table2array(post_DMTS_obj((post_DMTS_obj.Var4(:)==ffcc_array{1} | post_DMTS_obj.Var4(:)==ffcc_array{2}),2));
    onset_post_DMTS_city = table2array(post_DMTS_obj((post_DMTS_obj.Var4(:)==ffcc_array{3} | post_DMTS_obj.Var4(:)==ffcc_array{4}),2));
    
    %% 모든 condition들을 cell array로 묶기
    names = {'corr_match', 'corr_nonmatch', 'incorr_match', 'incorr_nonmatch', 'onset_pre_DMTS_forest', 'onset_pre_DMTS_city', 'onset_post_DMTS_forest', 'onset_post_DMTS_city'};
    onsets = {corr_match, corr_nonmatch, incorr_match, incorr_nonmatch, onset_pre_DMTS_forest, onset_pre_DMTS_city, onset_post_DMTS_forest, onset_post_DMTS_city};
    durations = {4, 4, 4, 4, 4, 4, 4, 4}; % 모든 duration은 4초
    
    multiple_conditions = {names, onsets, durations};

    if is_save_output == 1
%         save(fullfile(path_out{3}, [c_sbj, '_onset_pre_DMTS_forest.mat']),"onset_pre_DMTS_forest",'-mat')
%         save(fullfile(path_out{3}, [c_sbj, '_onset_pre_DMTS_city.mat']),"onset_pre_DMTS_city",'-mat')
%         save(fullfile(path_out{3}, [c_sbj, '_onset_post_DMTS_forest.mat']),"onset_post_DMTS_forest",'-mat')
%         save(fullfile(path_out{3}, [c_sbj, '_onset_post_DMTS_city.mat']),"onset_post_DMTS_city",'-mat')
        writetable(sbj_info_file,[path_out{2} '\sbj_info.xlsx']);
        save(fullfile(path_out{3},[c_sbj, '_multiple_conditions.mat']), 'names', 'onsets', 'durations');
    end
    
    %% Movement Regressor 저장
%% Path containing data
% path for confounding factors
conf_file = dir(fullfile(path{4}, c_sbj,'func', ['*task-', task_name, '*confounds*.tsv']));   % find the file

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
    if mod(sbj_i-1, 4) == 0
        if is_open_plot == 1
            f = figure('Position', [1500 500 1000 600]);
        else
            f = figure('Position', [1500 500 1000 600], 'Visible', 'off');
        end
    end

    % Select the subplot
    subplot(2, 2, mod(sbj_i-1, 4)+1);


    % Generate the plot for the current subject
        hold on
        colororder({'#0072BD','#000000'})

% RT plot : left
X=1:height(event_table);
yyaxis left
title([c_sbj ': RT & Correctness'],"FontSize",18,"FontWeight","bold")

% RT plot
plot(event_table.RT,'Color','#0072BD', 'LineWidth',1.8,'Marker','.','MarkerSize',20);
% Timeout threshold
yline(1.5,'-.','Timeout (>= 1.5s)','LabelHorizontalAlignment', 'center' ,'Color',"#0072BD",'LineWidth',1.2);

xlabel('Trial','FontSize',15,'FontWeight','bold')
ylabel('RT','FontSize',15,'FontWeight','bold')

yticks(0:0.2:1.8)
xlim([1 height(event_table)]);
ylim([-0.2 2]);
pbaspect([2 1 1]);


% correctness plot : right
yyaxis right
Y=event_table.Correct_Num';
Y(1,Y == 2) = 1;

stairs(Y);
coloringX = [X;X];
coloringY = [Y;Y];
CorPlot = area(coloringX([2:end end]),coloringY(1:end), 'FaceColor', 'k');
ylim([0 10]);

%legend
img = imread('correctness_legend.jpg');
image(img,'XData',[32 34],'YData',[3 0],'Clipping','off')

%timeout
coloringTO = [event_table.isTimeout'; event_table.isTimeout'];
area(coloringX([2:end end]), coloringTO(1:end), 'FaceColor',"#EDB120");
yticks([])


    % Save the plot for the current subject
    if is_save_output == 1
        saveas(gca, [path_out{4} '\individual\' c_sbj '_Performance.png']);
    end

    % Save the plot for every 4 subjects
%% 현재 켜져있는 figure 창 순서대로 저장하기 
% sub-01to04 이런식으로 4개씩 이름 지정

% Loop over all open figures by figure number
for i = 1:8
    % Get the figure handle
    f = figure(i);

    % Calculate the subject range for the current figure
    sbj_range = ((i-1)*4+1):i*4;

    % Create the figure name
    figName = ['sub-' sprintf('%02d', sbj_range(1)) 'to' sprintf('%02d', sbj_range(end)) '_Performance.png'];

    % Save the figure
    saveas(f, fullfile(path_out{4}, figName));
end

hold off
if is_open_plot == 0
    close all
end


end

    
    
    
    
