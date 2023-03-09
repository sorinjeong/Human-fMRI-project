LOG = log_data(1:300,:);
LogStr = string(LOG);

event_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "] " + optionalPattern(lettersPattern+": ")+optionalPattern(lettersPattern+":  ");

timeNevent = table; event_pre = []; time_pre = [];
for s =1:length(LogStr)
   [time,event] = get_time_event_from_log(LogStr(s));
   time_pre(s) = time;
   event_pre{s} = event;
end
processed_event = squeeze(split(string(event_pre), event_pat));
processed_event = [processed_event(:,2)];

timeNevent.time = time_pre';
timeNevent.event = processed_event;
