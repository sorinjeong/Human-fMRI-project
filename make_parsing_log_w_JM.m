LOG = log_data(1:300,:);
LogStr = string(LOG);

col_before_event = find(LogStr(250)==']',1,'last');
event = LogChar(col_before_event:end)
LogStr(250)

timeNevent = table
for s =1:length(LogStr)
[time,event] = get_time_event_from_log(LogStr(s));
if ~isempty(time) & ~isempty(event);
timeNevent.time(s) = time;
timeNevent.event{s} = event;
end
end

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

