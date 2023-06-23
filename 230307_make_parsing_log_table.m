LOG = log_data(1:300,:);


coordinate_pat = ["X=", "Y=", "Z="];
TF= contains(LOG,coordinate_pat);
find_true=find(TF);

first_T=min(find_true);
trimmed_log = LOG(first_T:end,1);
find_non_coordinate = trimmed_log(find(contains(trimmed_log,coordinate_pat)==0));

event_only_log = find_non_coordinate;

%%
timestamp = {};
br_timestamp_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "]";
for i = 1: length(str_event_only_log)
timestamp(i,1)=cellstr(extract(str_event_only_log(i),br_timestamp_pat));
end
str_timestamp = string(timestamp);
str_event_only_log = string(event_only_log);
timestamp_pat = digitsPattern(4) + "." + digitsPattern(2);

%%
message = ["Pause", 	"causeevent timeframe", 	"timespent",  "timespent_num", ...
    "Timeframe",     "relocation", 	"enter_arm",  "enter_arm_num",  ...
    "enter_gate_", "enter_gate_num",  "exit_gate_",  "exit_gate_num", ...
    "causeevent scanend", 	"task number",  "task_umber_num",  	...
    "start_control",  "start_control_num", 	"end_control",  "end_control_num", ...
    "correcttrials",  "correcttrials_num",  "answer",  "answer_num",  ...	
    "choice",  "choice_num",  "start_OCPR",  "start_OCPR_num", 	"end_OCPR", ...
    "end_OCPR_num",  "decision",  "decision_num", 	"period",  "period_num", ...	
    "Others",  "Note"];

        index_others = find(strcmp(message, 'Others'));
        index_note = find(strcmp(message, 'Note'));

var_oth = message(index_others);
var_not = message(index_note);
for j = 1:(length(message)-2)
    T =table;
    T.Properties.VariableNames = message(:)

    for i=1:length(str_event_only_log)
        logic = contains(str_event_only_log(i),message(j));

        message_N = char(strcat(message(j),'_num'));
        index_message_N = find(contains(message,message_N));

        if isempty(str_event_only_log(logic == 1)) == 0
            T.message(j) = cellstr(extract(str_timestamp(i),timestamp_pat));


            event_message = extractAfter(str_event_only_log(i),":  ");
            find_numeric = extract(event_message,digitsPattern);
            if isnumeric(find_numeric) == 1
                T.message(index_message_N) = find_numeric;
            end

        else 
            T.message(index_others)= cellstr(extract(str_timestamp(i), timestamp_pat));
            T.message(index_note) = cellstr(extractAfter(str_event_only_log(i),":  "));

        end

    end
%         var = [var; var(j)]
%         varo = [varo; var_oth]
%         varn = [varn; var_not]
% 
% 
% 
%     T = table(T(:,:), var);

end


% T = table(T(:,:),varo, varn);
% T.Properties.VariableNames = message(:)
T















