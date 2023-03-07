function [time,event] = get_time_event_from_log(LogStr)

brackets = find(LogStr=='[' | LogStr==']');
col_before_event = find(LogStr==':',1,'last');
LogChar = char(LogStr);

time = str2double(extractBetween(LogStr,"[","]"));
if ~isempty(time),time=time(1); end
event = LogChar(10:end);

end
