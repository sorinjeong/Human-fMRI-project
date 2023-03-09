function [time,event] = get_time_event_from_log(LogStr)

time = str2double(extractBetween(LogStr,"[","]"));
if ~isempty(time),time=time(1); end
event_log = LogStr;
event_pat = "["+digitsPattern(4) + "." + digitsPattern(2) + "] " + optionalPattern(lettersPattern+": ")+optionalPattern(lettersPattern+":  ");

if ~startsWith(event_log,"["+digitsPattern)
    event_log = "[0000.00] ";
    time = 0;
end
event = event_log;


end
