clear all; clc;
addpath(genpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV'));

%% set path
log_path_in = '../data/data_bhv_raw';
log_path_out = '../data/data_bhv_log_table'; 
plot_path_out = '../data/data_bhv_plot';

%% import the log file
flag_log = dir(fullfile(log_path_in,'*Behavior.csv'));
file_list = arrayfun(@(x) readtable(fullfile(x.folder, x.name)), flag_log, 'uni',0);
sbj_info_file = readtable('../../OCAT_DIR/data/data_fmri_bids/participants.tsv','FileType','text');
 
%% 
n_sbj = 31;

for sbj_i = 1: n_sbj
    c_sbj = num2str(sbj_i, '%02.f');

%% make output directory
   path_out = {};
   path_out{end+1} = fullfile(log_path_out,'individual',c_sbj);
   path_out{end+1} = fullfile(log_path_out,'total');
   path_out{end+1} = fullfile(log_path_out,'obj_showing_time');
   path_out{end+1} = fullfile(plot_path_out,'performance');
     
    if ~exist(log_path_out,"dir")
      mkdir(log_path_out);mkdir(plot_path_out);
      for i=1:length(path_out)
         mkdir(path_out{i});
      end
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
    writetable(event_TR,[path_out{1} '\' c_sbj '_event_TR.csv']);
    % remove MR row
    sbj_events(event_MR,:) = [];


    %% save PV task log
    PV_boundary = find(contains(sbj_events.Var1(:),'OCP'));
    event_pre_PV = sbj_events(PV_boundary(1):PV_boundary(2),:);
    event_post_PV = sbj_events(PV_boundary(3):PV_boundary(4),:);

    writetable(event_pre_PV,[path_out{1} '\' c_sbj '_event_pre_PV.csv']);
    writetable(event_post_PV,[path_out{1} '\' c_sbj '_event_post_PV.csv']);


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


type_log = find(contains(sbj_events.Var3(:),'Type'));
for i = 6:10
    event_struct.(var_name{i}) = num2str(sbj_events.Var4(type_log), i-5) - '0';
end


%% parsing events
event_name = string(sbj_events.Var1(:));
% Lap & Trial
time_lap_start = sbj_events(event_name=="LapStart",[2 4]);
for i=1:height(sbj_events)
    if event_name(i)=="TrialStart"
        event_struct.TrialStart(end+1) = sbj_events{i,2};
    elseif event_name(i)=="Trial"
        event_struct.Trial(end+1) = sbj_events{i,2};
        event_struct.Lap_Trial(end+1) = mod(length(event_struct.Trial), 4) + 1;
        lapidx = find(event_struct.TrialStart(end) > time_lap_start.Var2(1), 1, 'last');
        if ~isempty(lapidx); event_struct.Lap(end+1) = time_lap_start.Var4(lapidx);end

   % Choice
    else
        if contains(event_name(i), "Choice"+("A"|"B"))
        event_struct.Choice_txt{end+1} = string(extractAfter(event_name(i),"Choice"));
   % correctness, RT, isTimeout
        elseif event_name(i)=="Decision"
        %correct Timeout error
                if sbj_events{i,4} > 1.5 && sbj_events{i,2} ~= 2
                    sbj_events{i,2} = 2;
                end

             %correctness
             event_struct.Correct_Num(end+1) = sbj_events{i,2};
             %RT
             event_struct.RT(end+1) = sbj_events{i,4};

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
             event_struct.(event_name(i)){end+1} = sbj_events{i,2};
        end
    end

    %choice missing (choiceOn과 objOff 이벤트 개수와 choice 개수가 다를 경우 missing 삽입)
    if length(event_struct.ChoiceOn) == length(event_struct.ObjOff) && length(event_struct.ChoiceOn) ~= length(event_struct.Choice_txt)
    event_struct.Choice_txt{end+1} = {missing};
    end
end

% Choice txt to Number
event_struct.Choice_Num = NaN(height(event_struct.Choice_txt), 1);  % Initialize Choice_Num with NaN
event_struct.Choice_Num(strcmp(event_struct.Choice_txt, "A")) = 1;
event_struct.Choice_Num(strcmp(event_struct.Choice_txt, "B")) = 2;

% Context Num to txt
event_struct.Context_txt = replace(num2str(event_struct.Context_Num), {"1", "2"}, {"F", "C"});

%% Making a Table
event_struct = orderfields(event_struct,var_name);
event_table = struct2table(event_struct);

%% save the table
%individual
writetable(event_table,[path_out{1} '\' c_sbj '_event_table.csv']);
% total table
    %TR
    writetable(event_TR,[path_out{2} '\event_TR.csv'],'Sheet',c_sbj);
    %PV
    writetable(event_pre_PV,[path_out{2} '\event_pre_PV.csv'],'Sheet',c_sbj);
    writetable(event_post_PV,[path_out{2} '\event_post_PV.csv'],'Sheet',c_sbj);
    %events
    writetable(event_table,[path_out{2} '\event_table.csv'],'Sheet',c_sbj);


% correct_regressor; object가 켜진 시간, for GLM
corr = event_struct.ObjOn(event_struct.Correct_Num == 1);
incorr = event_struct.ObjOn(event_struct.Correct_Num ~= 1);

save([path_out{3} '\' c_sbj '_corr'],"corr",'-mat')
save([path_out{3} '\' c_sbj '_incorr'] ,"incorr",'-mat')



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
writetable(event_table,[path_out{1} '\' c_sbj '_event_numeric_table.csv']);
% save /all subjects
all_sbj_events = [all_sbj_events;event_numeric];


%% %%%%%%% performance plot (accuracy, RT) %%%%%%%%%
f = figure('Position', [1500 500 1000 600]);
hold on
colororder({'#0072BD','#000000'})

% RT plot : left
yyaxis left
title([c_sbj ': RT & Correctness'],"FontSize",18,"FontWeight","bold")

% RT plot
plot(event_table,"RT",'Color','#0072BD', 'LineWidth',1.8,'Marker','.','MarkerSize',20);
% Timeout threshold
yline(1.5,'-.','Timeout','LabelHorizontalAlignment', 'center' ,'Color',"#0072BD",'LineWidth',1.2);

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

area(repmat(1:height(event_table), 2, 1), repmat(Y, 2, 1), 'FaceColor', 'k');
ylim([0 10]);

%legend
img = imread('correctness_legend.jpg');
image(img,'XData',[32 34],'YData',[3 0],'Clipping','off')

%timeout
area([1:height(event_table);1:height(event_table)], [event_table.isTimeout;event_table.isTimeout], 'FaceColor',"#EDB120");
yticks([])

hold off
saveas(gcf,[path_out{4} '\' c_sbj '_Performance.png'])
close(f)


end
writetable(all_sbj_events,[path_out{2} '\all_sbj_events.csv']);




































