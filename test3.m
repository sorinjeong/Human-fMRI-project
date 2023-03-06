LOG = log_data(1:300,:);


pat = ["X=", "Y=", "Z="];
TF= contains(LOG,pat);
find_T=find(TF);

first_T=min(find_T);
trimmed_log = LOG(first_T:end,1);

%숫자 있는거 모두 열추가
message = ['Pause' 	'causeevent timeframe' 	'timespent'  'timespent_num' ...
    'Timeframe'     'relocation' 	'enter_arm'  'enter_arm_num'  ...
    'enter_gate_' 'enter_gate_num'  'exit_gate_'  'exit_gate_num' ...
    'causeevent scanend' 	'task number'  'task_umber_num'  	...
    'start_control'  'start_control_num' 	'end_control'  'end_control_num' ...
    'correcttrials'  'correcttrials_num'  'answer'  'answer_num'  ...	
    'choice'  'choice_num'  'start_OCPR'  'start_OCPR_num' 	'end_OCPR' ...
    'end_OCPR_num'  'decision'  'decision_num' 	'period'  'period_num' ...	
    'Others'  'Note'];

M=[];
time_sec = strcat("[","%f","]");
error_discription = strcat(":  ","%f");
% each_message = extractAfter(trimmed_log(i),":  ")
% find_string = sscanf(each_message,'%u')
str_trimmed_log = string(trimmed_log);

for i=1:length(trimmed_log)
    for j=1:length(message)
        logic = contains(trimmed_log(i),message(j));
        M(1,j) = message(1,j);
        N = strcat(message(1,j),'_num');

        if logic == 1
           M(end+1, j)= sscanf(str_trimmed_log(i), time_sec);
          
           each_message = extractAfter(str_trimmed_log(i),":  ");
           find_string = extractAfter(each_message,'%c');
           if isnumeric(find_string) == 1
               M(end+1,message(N)) = find_string;
           end

        else 
%              trimmed_log_char = cast(trimmed_log, 'char')
            M(end+1,message('Others')) = sscanf(str_trimmed_log(i), time_sec);
            M(end,message('Note')) = extractAfter(str_trimmed_log(i),":  ");
        end 
    end
end

Error_Table = M;
Error_Table