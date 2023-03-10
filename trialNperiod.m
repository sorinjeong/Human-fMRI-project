FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
    'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};

%%
for fi = 1:numel(FileList)
    filenameo=FileList{fi};
    filefolder= ['Y:\EPhysRawData\fmri_oppa_analysis\' filenameo];
    addpath(filefolder)

    %% Import the file
    Root =[ filefolder '\'];
    fileToRead1=strcat(filenameo,'B.log');
    DELIMITER = ' ';
    HEADERLINES = 50000;
    no_trials=80;

    rawData1 = importdata([Root fileToRead1], DELIMITER, HEADERLINES);
    log_data=rawData1;

    %% Log file to string
    LogStr = string(log_data(:));

EventName = ["Pause", "causeevent timeframe","timespent", "Timeframe", ...
    "relocation", "enter_arm", "enter_gate_", "exit_gate_", "causeevent scanend",...
    "task number", "start_control", "end_control", "correcttrials", "answer", ...	
    "choice", "start_OCPR", "end_OCPR", "decision", "period"];

VarName = ["Pause", "CTF", "timespent", "Timeframe", "reloc", "arm",...
    "entgate", "exgate", "CSC", "TN", "scont", "endcont", "correct",...
    "answer", "choice", "sOCPR", "eOCPR", "decision", "period"];

    %% make time and event table using function
    addpath 'C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange'
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

    %% test whether the log contains each event name
    coordinate_pat = ["X=", "Y=", "Z="];

    S=struct;
    Others=[];E=[];
for i = 1:height(timeNevent)
    if contains(timeNevent.event(i),coordinate_pat) ==1
        continue;
    else 
    for j = 1: length(EventName)
    contain_TF = contains(timeNevent.event{i},optionalPattern(digitsPattern) + EventName{j} + optionalPattern(digitsPattern), "IgnoreCase",true);

     if contain_TF == 1
        if ~isfield(S,VarName{j})
            S.(VarName{j}) = [];
        end 
         digit = [];
         digit = extract(timeNevent.event(i),digitsPattern);
       if ~isempty(digit)
           S.(VarName{j})= [S.(VarName{j}); timeNevent.time(i), extractBefore(timeNevent.event(i), EventName{j})+extractAfter(timeNevent.event(i), EventName{j})];   
       else
           S.(VarName{j})= [S.(VarName{j}); timeNevent.time(i), ""];
       end

       break;

      else 
         continue;
      end
           
   end

   if ~contain_TF
    Others = [Others; timeNevent.time(i), timeNevent.event(i)];
    S.Others = Others; 
   end

  end 
end 

S.eOCPR = [S.eOCPR; timeNevent.time(end), timeNevent.event(end)];
save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo], "S")



%% make trial N period

Str=struct('Time',[],'Trial',[],'Period',[],'Event',[]);

var_start = ["scont", "sOCPR"];
var_end = ["endcont", "eOCPR"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% e=1 control only | e=2 OCPR only | e=3 BOTH

for c=1:length(S.scont)
%     e=1; % control only
%     e=2; % OCPR only
    for e=1:length(var_start) % both
        A=S.(var_start{e}); B=S.(var_end{e});

    startTime = A(c,1);
    endTime = B(c,1);
    timeBetween_p = startTime <= S.period(1:end,1) & S.period(1:end,1) <= endTime;
    timeBetween_d = startTime <= S.decision(1:end,1) & S.decision(1:end,1) <= endTime;
    [tp_r, tp_c] = find(timeBetween_p);
    [td_r, td_c] = find(timeBetween_d);

%% start
    S_Time = startTime;
    S_Trial = A(c,2);
    S_Period = "1";
    S_Event = "start";

%% make struct

if isempty(Str)
    Str.Time= S_Time;
    Str.Trial= S_Trial;
    Str.Period= S_Period;
    Str.Event= S_Event;
else

    Str.Time= [Str.Time; S_Time];
    Str.Trial= [Str.Trial; S_Trial];
    Str.Period= [Str.Period; S_Period];
    Str.Event= [Str.Event; S_Event];
end
%% dicision
% if ~isempty(td_r)    
    D_Time = unique(S.decision(td_r,1));
    D_Trial = Str.Trial(end);
    D_Period = Str.Period(end);
    D_Event = "decision";



    %% period
    P_Time = S.period(tp_r,1);
    P_Trial = Str.Trial(end);
    P_Period = extract(S.period(tp_r,2), digitsPattern);
    P_Event = "start";



%% dicision vs period
if ~isempty(D_Time)

for i=1:length(P_Time)
    for j=1:length(D_Time)
        if P_Time(i) > D_Time(j)
            Str.Time= [Str.Time; D_Time(j)];
%             D_Time(j) = NaN;
            Str.Trial = [Str.Trial; D_Trial];
                        S.decision(td_r,3) = D_Trial;
            Str.Period = [Str.Period; D_Period];
                        S.decision(td_r,4) = D_Period;
            Str.Event = [Str.Event; D_Event];
        else 
            P_Time(i) < D_Time(j);
            break;
        end
    end
    Str.Time = [Str.Time; P_Time(i)];
            Str.Trial = [Str.Trial; P_Trial];
            Str.Period = [Str.Period; P_Period(i)];
            Str.Event = [Str.Event; P_Event];
end

if length(P_Time) < length(D_Time)
    for a=length(P_Time)+1:length(D_Time)
        Str.Time=[Str.Time; D_Time(a)];
            Str.Trial = [Str.Trial; D_Trial];
                        S.decision(td_r,3) = D_Trial;
            Str.Period = [Str.Period; D_Period ];
                        S.decision(td_r,4) = D_Period;
            Str.Event = [Str.Event; D_Event];
    end
elseif length(P_Time) == length(D_Time)
    Str.Time=[Str.Time; D_Time(length(D_Time))];
            Str.Trial = [Str.Trial; D_Trial];
                         S.decision(td_r,3) = D_Trial;
            Str.Period = [Str.Period; D_Period ];
                         S.decision(td_r,4) = D_Period;
            Str.Event = [Str.Event; D_Event];
end
end
    %% end
    E_Time = endTime;
    E_Trial = Str.Trial(end);
    E_Period = Str.Period(end);
    E_Event = "end";

    Str.Time= [Str.Time; E_Time];
    Str.Trial= [Str.Trial; E_Trial];
    Str.Period= [Str.Period; E_Period];
    Str.Event= [Str.Event; E_Event];

    end  % both
end
%% make a table
trialNperiod = table(Str.Time, Str.Trial, Str.Period, Str.Event);
trialNperiod.Properties.VariableNames = ["Time","Trial","Period","Event"];

% tablename = 'trialNperiod_control_only'; % control only
% tablename = 'trialNperiod_OCPR_only'; % OCPR only
tablename = 'trialNperiod_BOTH'; % both

save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo '_' tablename], "trialNperiod");

end



