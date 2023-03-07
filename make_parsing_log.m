LOG = log_data(1:300,:);


coordinate_pat = ["X=", "Y=", "Z="];
TF= contains(LOG,coordinate_pat);
find_T=find(TF);

first_T=min(find_T);
trimmed_log = LOG(first_T:end,1);

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

%character, cell
message = {'Pause', 	'causeevent timeframe', 	'timespent',  'timespent_num', ...
    'Timeframe',     'relocation', 	'enter_arm',  'enter_arm_num',  ...
    'enter_gate_', 'enter_gate_num',  'exit_gate_',  'exit_gate_num', ...
    'causeevent scanend', 	'task number',  'task_umber_num',  	...
    'start_control',  'start_control_num', 	'end_control',  'end_control_num', ...
    'correcttrials',  'correcttrials_num',  'answer',  'answer_num',  ...	
    'choice',  'choice_num',  'start_OCPR',  'start_OCPR_num', 	'end_OCPR', ...
    'end_OCPR_num',  'decision',  'decision_num', 	'period',  'period_num', ...	
    'Others',  'Note'}



M={};
M(1,:)=message(:)
% time_sec = strcat("[","%f","]");
% error_discription = strcat(":  ","%f");
timestamp_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "]";
% each_message = extractAfter(trimmed_log(i),":  ")
% find_numeric = sscanf(each_message,'%u')
str_trimmed_log = string(trimmed_log);

for i=1:length(str_trimmed_log)
    for j=1:length(message)
        logic = contains(str_trimmed_log(i),message(j));
        message_N = strcat(message(1,j),'_num');
        index_message_N = find(strcmp(message, message_N));

        if logic == 1
%            M(end+1, j)= sscanf(str_trimmed_log(i), time_sec);
            M(end+1, j) = cellstr(extract(str_trimmed_log(i),timestamp_pat));
          
              event_message = extractAfter(str_trimmed_log(i),":  ");
%              find_numeric = extractAfter(event_message,'%c');
              num_pat = digitBoundary("end");
              find_numeric = extract(event_message,num_pat);

           if isnumeric(find_numeric) == 1
               M(end,index_message_N) = cellstr(find_numeric);
           end

        else 
%              trimmed_log_char = cast(trimmed_log, 'char')
%             M(end+1,message('Others')) = sscanf(str_trimmed_log(i), time_sec);
            index_others = find(strcmp(message, 'Others'));
            index_note = find(strcmp(message, 'Note'));
            M(end+1,index_others) = cellstr(extract(str_trimmed_log(i), timestamp_pat));
            M(end,index_note) = cellstr(extractAfter(str_trimmed_log(i),":  "));
        end 
    end
end

br_pat = ["[", "]"];
M_final = erase(M{:},br_pat);

Event_Table = M_final;
Event_Table
