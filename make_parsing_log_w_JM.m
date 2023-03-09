LOG = log_data(1:300,:);
LogStr = string(LOG);

event_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "] " + optionalPattern(lettersPattern+": ")+optionalPattern(lettersPattern+":  ");

timeNevent = table;
event_pre = [];
time_pre = [];
% for s = 1:length(LogStr)
%     [time,event] = get_time_event_from_log(LogStr(s));
        for t =1:length(LogStr)
         [time,event] = get_time_event_from_log(LogStr(t));
         time_pre(t) = time;
         event_pre{t} = event;
        end
         processed_event = squeeze(split(string(event_pre), event_pat));
         processed_event = [processed_event(:,2)];
     

%     if ~isempty(time)
        timeNevent.time = time_pre';
%     end
    timeNevent.event = processed_event;
% end
%%


var= [];T=table;
event_name = ["Pause", "causeevent_timeframe"];
    for j = 1: length(event_name)
        for i = 1:height(timeNevent)
            event_col = timeNevent{i,'event'};
            contain_TF = strcmp(event_col,"*"+event_name(j)+"*");

        if contain_TF == 0;
           var = timeNevent{i,"time"} ;


        end
        end

var(j) = [var(j); var];
tablevariables = [];
tablevariables(1,j) = var(j)
T = table(VariableNames);
    T.Properties.VariableNames = event_name(i);



end

