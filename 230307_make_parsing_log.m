LOG = log_data(1:300,:);


coordinate_pat = ["X=", "Y=", "Z="];
TF= contains(LOG,coordinate_pat);
find_true=find(TF);

first_T=min(find_true);
trimmed_log = LOG(first_T:end,1);
% find_coordinate = trimmed_log(find(contains(trimmed_log,coordinate_pat)));
find_non_coordinate = trimmed_log(find(contains(trimmed_log,coordinate_pat)==0));

event_only_log = find_non_coordinate;
%%

%숫자 있는거 모두 열추가
%string
% message = ["Pause", 	"causeevent timeframe", 	"timespent",  "timespent_num", ...
%     "Timeframe",     "relocation", 	"enter_arm",  "enter_arm_num",  ...
%     "enter_gate_", "enter_gate_num",  "exit_gate_",  "exit_gate_num", ...
%     "causeevent scanend", 	"task number",  "task_umber_num",  	...
%     "start_control",  "start_control_num", 	"end_control",  "end_control_num", ...
%     "correcttrials",  "correcttrials_num",  "answer",  "answer_num",  ...	
%     "choice",  "choice_num",  "start_OCPR",  "start_OCPR_num", 	"end_OCPR", ...
%     "end_OCPR_num",  "decision",  "decision_num", 	"period",  "period_num", ...	
%     "Others",  "Note"];

% %character, cell
% message = {'Pause', 	'causeevent timeframe', 	'timespent',  'timespent_num', ...
%     'Timeframe',     'relocation', 	'enter_arm',  'enter_arm_num',  ...
%     'enter_gate_', 'enter_gate_num',  'exit_gate_',  'exit_gate_num', ...
%     'causeevent scanend', 	'task number',  'task_umber_num',  	...
%     'start_control',  'start_control_num', 	'end_control',  'end_control_num', ...
%     'correcttrials',  'correcttrials_num',  'answer',  'answer_num',  ...	
%     'choice',  'choice_num',  'start_OCPR',  'start_OCPR_num', 	'end_OCPR', ...
%     'end_OCPR_num',  'decision',  'decision_num', 	'period',  'period_num', ...	
%     'Others',  'Note'};

%cell_ver2
message = strsplit('Pause,causeevent timeframe,timespent,timespent_num,Timeframe,relocation,enter_arm,enter_arm_num,enter_gate_,enter_gate_num,exit_gate_,exit_gate_num,causeevent scanend,task number,task_umber_num,start_control,start_control_num,end_control,end_control_num,correcttrials,correcttrials_num,answer,answer_num,choice,choice_num,start_OCPR,start_OCPR_num,end_OCPR,end_OCPR_num,decision,decision_num,period,period_num,Others,Note',',');


timestamp = {};
br_timestamp_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "]";
for i = 1: length(str_event_only_log)
timestamp(i,1)=cellstr(extract(str_event_only_log(i),br_timestamp_pat));
end
str_timestamp = string(timestamp);

M={};
M(1,:)=message(:);
% time_sec = strcat("[","%f","]");
% error_discription = strcat(":  ","%f");
% each_message = extractAfter(event_only_log(i),":  ")
% find_numeric = sscanf(each_message,'%u')
str_event_only_log = string(event_only_log);



%% 


for i=1:length(str_event_only_log)
    for j=1:length(message)
        logic = contains(str_event_only_log(i),message(j));

        timestamp_pat = digitsPattern(4) + "." + digitsPattern(2);
        message_N = char(strcat(message(1,j),'_num'));
        index_message_N = find(contains(message,message_N));
        index_others = find(strcmp(message, 'Others'));
        index_note = find(strcmp(message, 'Note'));
       

        if isempty(str_event_only_log(logic == 1)) == 0
%            M(end+1, j)= sscanf(str_trimmed_log(i), time_sec);
            M(i+1, j) = cellstr(extract(str_timestamp(i),timestamp_pat));

        else 
%             trimmed_log_char = cast(trimmed_log, 'char')
%             M(end+1,message('Others')) = sscanf(str_trimmed_log(i), time_sec);

            M(i+1,index_others) = cellstr(extract(str_timestamp(i), timestamp_pat));
            M(i+1,index_note) = cellstr(extractAfter(str_event_only_log(i),":  "));
        end
        
          %%
              event_message = extractAfter(str_event_only_log(i),":  ");
%              find_numeric = extractAfter(event_message,'%c');
              find_numeric = extract(event_message,digitsPattern);

           if isnumeric(find_numeric) == 1
               M(i+1,index_message_N) = find_numeric;
           end

       
    end
end


Event_Table = M;
Event_Table
