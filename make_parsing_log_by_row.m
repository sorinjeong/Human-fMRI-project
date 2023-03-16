% FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
%     'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};

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
    addpath 'C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange'
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
mkdir(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed' filenameo]); % OCPR / control / BOTH 중 처음 한번만 run!
save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo], "S");



%% make trial N period

Str=struct('Time',[],'Trial',[],'Period',[],'Event',[],'Note',[]);

var_start = ["scont", "sOCPR"];
var_end = ["endcont", "eOCPR"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% e=1 control only | e=2 OCPR only | e=3 BOTH
c=missing;e=missing;
for c=1:length(S.scont)
%     e=1; % control only
%     e=2; % OCPR only
    for e=1:length(var_start) % both
        A=S.(var_start{e}); B=S.(var_end{e}); [o,oc] = find(S.Others == "9");

if c ~= length(A)
    startLog = str2double(A(c,3));
    nextstartLog = str2double(A(c+1,3));
    logBetween_p = startLog <= str2double(S.period(:,3)) & str2double(S.period(:,3)) < nextstartLog;
    logBetween_d = startLog <= str2double(S.decision(:,3)) & str2double(S.decision(:,3)) < nextstartLog;
    logBetween_cd = startLog <= str2double(S.caus_decision(:,3)) & str2double(S.caus_decision(:,3)) < nextstartLog;
    [tp_r, tp_c] = find(logBetween_p);
    [td_r, td_c] = find(logBetween_d);
    [tcd_r, tcd_c] = find(logBetween_cd);

    logBetween_c = startLog <= str2double(S.choice(:,3)) & str2double(S.choice(:,3)) < nextstartLog;
    logBetween_a = startLog <= str2double(S.answer(:,3)) & str2double(S.answer(:,3)) < nextstartLog;
    logBetween_o = startLog <= str2double(S.Others(o,3)) & str2double(S.Others(o,3)) < nextstartLog;

    [lc_r, lc_c] = find(logBetween_c);
    [la_r, la_c] = find(logBetween_a);
    [lo_r, lo_c] = find(logBetween_o);

elseif c == length(A)
    startLog = str2double(A(c,3));
    logAfter_p = startLog <= str2double(S.period(:,3));
    logAfter_d = startLog <= str2double(S.decision(:,3));
    logAfter_cd = startLog <= str2double(S.caus_decision(:,3));

    [tp_r, tp_c] = find(logAfter_p);
    [td_r, td_c] = find(logAfter_d);
    [tcd_r, tcd_c] = find(logAfter_cd);

    logAfter_c = startLog <= str2double(S.choice(:,3));
    logAfter_a = startLog <= str2double(S.answer(:,3));
    logAfter_o = startLog <= str2double(S.Others(o,3));

    [lc_r, lc_c] = find(logAfter_c);
    [la_r, la_c] = find(logAfter_a);
    [lo_r, lo_c] = find(logAfter_o);
end

%% start
    S_Time = startLog;
    S_Trial = A(c,2);
    S_Period = "1";
    S_Event = "start";
    S_Note = missing;

%% make struct

if isempty(Str)
    Str.Time= S_Time;
    Str.Trial= S_Trial;
    Str.Period= S_Period;
    Str.Event= S_Event;
    Str.Note = S_Note;
else

    Str.Time= [Str.Time; S_Time];
    Str.Trial= [Str.Trial; S_Trial];
    Str.Period= [Str.Period; S_Period];
    Str.Event= [Str.Event; S_Event];
    Str.Note= [Str.Note; S_Note];
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
    E_Time = nextstartLog;
    E_Trial = S_Trial;
    E_Event = "end";
    E_Log = str2double(B(c,3));

%% correctness answer and choice #1
C_Log=[]; A_Log=[];

    C_Time = str2double(S.choice(lc_r,1)');
    C_Note = str2double(S.choice(lc_r,2)');
    A_Time = str2double(S.answer(la_r,1)');
    A_Note = str2double(S.answer(la_r,2)');
    O_Time = str2double(S.Others(lo_r,1)');
    O_Note = str2double(S.Others(lo_r,2)');


for i=1:height(lc_r)
    for j=1:height(la_r)
        if floor(C_Time(i)) == floor(A_Time(j))
            C_Log = [C_Log, str2double(S.choice(lc_r(i),3))];
            A_Log = [A_Log, str2double(S.answer(la_r(j),3))];
            break;
        else 
            for k=1:length(lo_r)
                if floor(C_Time(i)) == floor(O_Time(k))
                    C_Log = [C_Log, str2double(S.choice(lc_r(i),3))];
                    A_Log = [A_Log, str2double(S.Others(lo_r(k),3))];  
                    break;
                end
            end
        end
        if ~isempty(C_Log)
        if C_Log(end) == str2double(S.choice(lc_r(i),3))
            break;
        end
        end
    end
end


%% new strategy decision vs period

log_order = sort([P_Log D_Log CD_Log E_Log A_Log]);

for i=1:length(log_order)

    P_index = find(log_order(i) == P_Log);
    D_index = find(log_order(i) == D_Log);
    CD_index = find(log_order(i) == CD_Log);
    E_index = find(log_order(i) == E_Log);
    A_index = find(log_order(i) == A_Log);

    if ~isempty(P_index)
        Str.Time = [Str.Time; P_Time(P_index)];
        Str.Trial = [Str.Trial; P_Trial];
        Str.Period = [Str.Period; P_Period(P_index)];
        Str.Event = [Str.Event; P_Event];
        Str.Note = [Str.Note; missing];

    elseif ~isempty(D_index)
        Str.Time= [Str.Time; D_Time(D_index)];
        Str.Trial = [Str.Trial; missing];
        Str.Period = [Str.Period; missing];
        Str.Event = [Str.Event; D_Event];
        Str.Note = [Str.Note; missing];

    elseif ~isempty(CD_index)
        Str.Time= [Str.Time; CD_Time(CD_index)];
        Str.Trial = [Str.Trial; missing];
        Str.Period = [Str.Period; missing];
        Str.Event = [Str.Event; CD_Event];
        Str.Note = [Str.Note; missing];

    elseif ~isempty(E_index)
        Str.Time= [Str.Time; E_Time];
        Str.Trial= [Str.Trial; missing];
        Str.Period= [Str.Period; missing];
        Str.Event= [Str.Event; E_Event];
        Str.Note = [Str.Note; missing];

    elseif ~isempty(A_index)
        answer_log = str2double(S.answer(la_r,3));
        others_log = str2double(S.Others(lo_r,3));
        if ~isempty(A_Log(A_index) == answer_log)
          Str.Time= [Str.Time; A_Time(find(la_r(A_Log(A_index) == answer_log)))];
          Str.Event = [Str.Event; "answer"];
          Str.Note = [Str.Note; A_Note(find(la_r(A_Log(A_index) == answer_log)))];
        elseif ~isempty(A_Log(A_index) == others_log)
          Str.Time= [Str.Time; O_Time(find(lo_r(A_Log(A_index) == others_log)))];
          Str.Event = [Str.Event; "miss_answer9"];
          Str.Note = [Str.Note; O_Note(find(lo_r(A_Log(A_index) == others_log)))];
        end
          Str.Time= [Str.Time; C_Time(floor(C_Time)==floor(str2double(Str.Time(end))))];
          Str.Event = [Str.Event; "choice"];
          Str.Note = [Str.Note; C_Note(floor(C_Time)==floor(str2double(Str.Time(end))))];
          Str.Trial = [Str.Trial; missing; missing];
          Str.Period = [Str.Period; missing; missing];
    end 
end
Str.Trial = fillmissing(Str.Trial,'previous');
Str.Period = fillmissing(Str.Period,'previous');
    end  % both
end

%% make a table
AnswerChoice = table(str2double(Str.Time), str2double(Str.Trial), str2double(Str.Period), Str.Event, Str.Note);
AnswerChoice.Properties.VariableNames = ["Time","Trial","Period","Event", "Note_direction"];

% tablename = 'trialNperiod_control_only'; % control only
% tablename = 'trialNperiod_OCPR_only'; % OCPR only
tablename = 'AnswerChoice'; % both

save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo '_' tablename], "AnswerChoice");


%%%%%%%%%%%%%%%%% both일 경우에만 아래 table 생성!%%%%%%%%%%%%%%%%%%%%%
%% S.decision + [trial period]

Str.decision = S.decision;
for i=1:height(AnswerChoice)
   for j=1:length(Str.decision)
     if AnswerChoice.Time(i) == str2double(Str.decision(j,1))
        Str.decision(j,4) = AnswerChoice.Trial(i);
        Str.decision(j,5) = AnswerChoice.Period(i);
     end
  end
end
% fillmissing(Str.decision,'constant',NaN)

Decision_table = splitvars(table(str2double(Str.decision)));
Decision_table.Properties.VariableNames = ["Time","ordi_num","Log","Trial","Period"];
save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo '_Decision'], "Decision_table");

%% S.caus_decision + [trial period]

Str.caus_decision = S.caus_decision;
for i=1:height(AnswerChoice)
   for j=1:length(Str.caus_decision)
      if AnswerChoice.Time(i) == str2double(Str.caus_decision(j,1))
        Str.caus_decision(j,4) = AnswerChoice.Trial(i);
        Str.caus_decision(j,5) = AnswerChoice.Period(i);
      end
   end
end

CausDecision_table = splitvars(table(str2double(Str.caus_decision)));
CausDecision_table.Properties.VariableNames = ["Time","ordi_num","Log","Trial","Period"];
CausDecision_table = removevars(CausDecision_table,"ordi_num");
save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo '_CausDecision'], "CausDecision_table");

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
save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo '_Pause'], "Pause_table");



end


