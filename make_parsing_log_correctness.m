FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
    'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};

FileList = {'CL130227_1'};
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
    "choice", "start_OCPR", "end_OCPR", "causeevent decision", "decision", "period"];

VarName = ["Pause", "CTF", "timespent", "Timeframe", "reloc", "arm",...
    "entgate", "exgate", "CSC", "TN", "scont", "endcont", "correct",...
    "answer", "choice", "sOCPR", "eOCPR", "caus_decision", "decision", "period"];

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
    Others=[];E=[];decision=[];
for i = 1:height(timeNevent)
    if contains(timeNevent.event(i),coordinate_pat) ==1 & ~contains(timeNevent.event{i}, "2nd decision")
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
                  S.(VarName{j})= [S.(VarName{j}); timeNevent.time(i), extractBefore(timeNevent.event(i), EventName{j})+extractAfter(timeNevent.event(i), EventName{j}), i];   
                 else
                   S.(VarName{j})= [S.(VarName{j}); timeNevent.time(i), "", i];
                end

              break;
            else 
             continue;
            end
           
        end

        if ~contain_TF
          Others = [Others; timeNevent.time(i), timeNevent.event(i),i];
          S.Others = Others; 
        end
    end 
end 

for i=1:height(S.period)
 S.period(i,2) = str2double(extractBefore(S.period(i,2),lettersPattern));
end
for i=1:height(S.decision)
 S.decision(i,2) = str2double(extractBefore(S.decision(i,2),lettersPattern));
end
for i=1:height(S.caus_decision)
 S.caus_decision(i,2) = str2double(extractBefore(S.caus_decision(i,2),lettersPattern));
end

S.eOCPR = [S.eOCPR; timeNevent.time(end), timeNevent.event(end),i];

%% save log data
% mkdir(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo]); % OCPR / control / BOTH 중 처음 한번만 run!
save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo '\' filenameo], "S");



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

    startTime = str2double(A(c,1));
    endTime = str2double(B(c,1));
    timeBetween_p = startTime <= str2double(S.period(1:end,1)) & str2double(S.period(1:end,1)) <= endTime;
    timeBetween_d = startTime <= str2double(S.decision(1:end,1)) & str2double(S.decision(1:end,1)) <= endTime;
    timeBetween_cd = startTime <= str2double(S.caus_decision(1:end,1)) & str2double(S.caus_decision(1:end,1)) <= endTime;
    [tp_r, tp_c] = find(timeBetween_p);
    [td_r, td_c] = find(timeBetween_d);
    [tcd_r, tcd_c] = find(timeBetween_cd);


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

%% period
    P_Time = S.period(tp_r,1)';
    P_Trial = S_Trial;
    P_Period = extract(S.period(tp_r,2), digitsPattern)';
    P_Event = "start";
    P_Log = str2double(S.period(tp_r,3)');

%% decision
% if ~isempty(td_r)    
    D_Time = unique(S.decision(td_r,1))';
    D_Event = "decision";
    D_Log=[];
    %% make D_Log
    repeat_log=[];
    for i=1:length(td_r)
    if ismissing(S.decision(td_r(i),2))
        repeat_log = [repeat_log, S.decision(td_r(i),3)];
        continue;
    else
        D_Log = [D_Log, S.decision(td_r(i),3)];
    end
    end
if ~isempty(repeat_log)
D_Log = [D_Log, repeat_log(1)];
end
D_Log = str2double(sort(D_Log));


%% causeevent decision

    CD_Time = unique(S.caus_decision(tcd_r,1))';
    CD_Event = "caus_decision";
    CD_Log = str2double(S.caus_decision(tcd_r,3)');

%% end
    E_Time = endTime;
    E_Trial = S_Trial;
    E_Event = "end";
    E_Log = str2double(B(c,3));

%% new strategy decision vs period

log_order = sort([P_Log D_Log CD_Log E_Log]);

for i=1:length(log_order)

    P_index = find(log_order(i) == P_Log);
    D_index = find(log_order(i) == D_Log);
    CD_index = find(log_order(i) == CD_Log);
    E_index = find(log_order(i) == E_Log);

    if ~isempty(P_index)
        Str.Time = [Str.Time; P_Time(P_index)];
        Str.Trial = [Str.Trial; P_Trial];
        Str.Period = [Str.Period; P_Period(P_index)];
        Str.Event = [Str.Event; P_Event];

    elseif ~isempty(D_index)
        Str.Time= [Str.Time; D_Time(D_index)];
        Str.Trial = [Str.Trial; Str.Trial(end)];
        Str.Period = [Str.Period; Str.Period(end)];
        Str.Event = [Str.Event; D_Event];

    elseif ~isempty(CD_index)
        Str.Time= [Str.Time; CD_Time(CD_index)];
        Str.Trial = [Str.Trial; Str.Trial(end)];
        Str.Period = [Str.Period; Str.Period(end)];
        Str.Event = [Str.Event; CD_Event];

    elseif ~isempty(E_index)
        Str.Time= [Str.Time; E_Time];
        Str.Trial= [Str.Trial; E_Trial];
        Str.Period= [Str.Period; Str.Period(end)];
        Str.Event= [Str.Event; E_Event];
    end 
end

    end  % both
end

%% make a table
trialNperiod = table(str2double(Str.Time), str2double(Str.Trial), str2double(Str.Period), Str.Event);
trialNperiod.Properties.VariableNames = ["Time","Trial","Period","Event"];

% tablename = 'trialNperiod_control_only'; % control only
% tablename = 'trialNperiod_OCPR_only'; % OCPR only
tablename = 'trialNperiod_BOTH'; % both

save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo '\' filenameo '_' tablename], "trialNperiod");


%%%%%%%%%%%%%%%%% both일 경우에만 아래 table 생성!%%%%%%%%%%%%%%%%%%%%%
%% S.decision + [trial period]

Str.decision = S.decision;
for i=1:height(trialNperiod)
   for j=1:length(Str.decision)
     if trialNperiod.Time(i) == str2double(Str.decision(j,1))
        Str.decision(j,4) = trialNperiod.Trial(i);
        Str.decision(j,5) = trialNperiod.Period(i);
     end
  end
end
% fillmissing(Str.decision,'constant',NaN)

Decision_table = splitvars(table(str2double(Str.decision)));
Decision_table.Properties.VariableNames = ["Time","ordi_num","Log","Trial","Period"];
save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo '\' filenameo '_Decision'], "Decision_table");

%% S.caus_decision + [trial period]

Str.caus_decision = S.caus_decision;
for i=1:height(trialNperiod)
   for j=1:length(Str.caus_decision)
      if trialNperiod.Time(i) == str2double(Str.caus_decision(j,1))
        Str.caus_decision(j,4) = trialNperiod.Trial(i);
        Str.caus_decision(j,5) = trialNperiod.Period(i);
      end
   end
end

CausDecision_table = splitvars(table(str2double(Str.caus_decision)));
CausDecision_table.Properties.VariableNames = ["Time","ordi_num","Log","Trial","Period"];
CausDecision_table = removevars(CausDecision_table,"ordi_num");
save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo '\' filenameo '_CausDecision'], "CausDecision_table");

%% S.Pause_table
Str.Pause = S.Pause;
for i = 1:2:height(S.Pause)
        Str.Pause(i,2) = "start";
end
for i = 2:2:height(S.Pause)
        Str.Pause(i,2) = "end";
end
Str.Pause(end,2) = "termination";
Pause_table = splitvars(table(Str.Pause));
Pause_table.Properties.VariableNames = ["Time","start/end","Log"];
save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo '\' filenameo '_Pause'], "Pause_table");

%% correctness answer and choice

[o,c] = find(S.Others == "9");

correctStruct = struct('Time',trialNperiod.Time,'Trial',trialNperiod.Trial,'Period',trialNperiod.Period,'Event',trialNperiod.Event,'Note',[NaN(height(trialNperiod),1)]);

for k=1:height(S.answer)
    for j=1:height(S.choice)
%         for a=1:length(o)
            for i=1:(height(correctStruct)+(height(S.choice)*2))
                if correctStruct.Time(i) == str2double(S.choice(j,1))
                    continue;
                elseif correctStruct.Time(i) < str2double(S.choice(j,1)) & correctStruct.Time(i+1) > str2double(S.choice(j,1))
                    for a=1:length(o)
                    if str2double(S.answer(k,1)) == str2double(S.choice(j,1))
                        correctStruct.Time = [correctStruct.Time(1:i); str2double(S.answer(k,1)); str2double(S.choice(j,1)); correctStruct.Time(i+1:end)];
                        correctStruct.Trial = [correctStruct.Trial(1:i); correctStruct.Trial(i); correctStruct.Trial(i); correctStruct.Trial(i+1:end)];
                        correctStruct.Period = [correctStruct.Period(1:i); correctStruct.Period(i); correctStruct.Period(i); correctStruct.Period(i+1:end)];
                        correctStruct.Event =[correctStruct.Event(1:i); "answer"; "choice"; correctStruct.Event(i+1:end)];
                        correctStruct.Note = [correctStruct.Note(1:i); str2double(S.answer(k,2)); str2double(S.choice(j,2)); correctStruct.Note(i+1:end)];
                        break;
                    elseif str2double(S.Others(o(a),1)) == str2double(S.choice(j,1))
                        correctStruct.Time = [correctStruct.Time(1:i); str2double(S.Others(o(a),1)); str2double(S.choice(j,1)); correctStruct.Time(i+1:end)];
                        correctStruct.Trial = [correctStruct.Trial(1:i); correctStruct.Trial(i); correctStruct.Trial(i); correctStruct.Trial(i+1:end)];
                        correctStruct.Period = [correctStruct.Period(1:i); correctStruct.Period(i); correctStruct.Period(i); correctStruct.Period(i+1:end)];
                        correctStruct.Event =[correctStruct.Event(1:i); "miss_answer9"; "choice"; correctStruct.Event(i+1:end)];
                        correctStruct.Note = [correctStruct.Note(1:i); 9; str2double(S.choice(j,2)); correctStruct.Note(i+1:end)];
end
                    end
%                 end
            end
        end
    end
end


correctTable = table(correctStruct.Time, correctStruct.Trial, correctStruct.Period, correctStruct.Event, correctStruct.Note);
correctTable.Properties.VariableNames = ["Time","Trial","Period","Event", "Note"];




end


