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
 S.period(i,2) = (extractBefore(S.period(i,2),lettersPattern));
end
for i=1:height(S.decision)
 S.decision(i,2) = (extractBefore(S.decision(i,2),lettersPattern));
end
for i=1:height(S.caus_decision)
 S.caus_decision(i,2) = (extractBefore(S.caus_decision(i,2),lettersPattern));
end

S.eOCPR = [S.eOCPR; timeNevent.time(end), timeNevent.event(end),i];

%% save log data
% mkdir(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo]); % OCPR / control / BOTH 중 처음 한번만 run!
save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo], "S");


%% arrange trial, event by log order

Str=struct('Time',[],'Trial',[],'Period',[],'Event',[],'Note',[]);
%% both: USE 'startNend'!!
startNend=[];

startNend = [(S.scont), (S.endcont); (S.sOCPR), (S.eOCPR)];
    startNend = sortrows(str2double(startNend),[3 6 1 4]);

%% make answer + others(answer_9) string
[o,oc] = find(S.Others == "9");
o=o(10<o);
answer_9_error = str2double([S.Others(o,:)]);
real_answer=sortrows([str2double(S.answer); answer_9_error],[3 1]);

Answers = string([real_answer(:,1), NaN(height(real_answer),3),real_answer(:,2:3)]);
for i=1:height(Answers)    
if str2double(Answers(i,1)) ~= answer_9_error(:,1)
    Answers(i,4) = "answers";
end
end
Answers(:,4) = fillmissing(Answers(:,4),'constant', "miss_answers9");

%% modify decision
cutDecision=S.decision;
    for i=1:height(S.decision)
        if ismissing(S.decision(i,2)) & ismissing(S.decision(i+1,2))
            cutDecision(i,:)=[];
        end
    end

    %% assemble total logs

% modify logs to 5column
modst = [S.scont(:,1:2); S.sOCPR(:,1:2)];
modst(:,3) = 1 ; modst(:,4) = "start"; modst(:,5) = missing; modst(:,6) = [S.scont(:,3); S.sOCPR(:,3)];
modend = [S.endcont(:,1:2); S.eOCPR(:,1:2)];
modend(:,3) = missing ; modend(:,4) = "end"; modend(:,5) = missing; modend(:,6) = [S.endcont(:,3); S.eOCPR(:,3)];
modper = [S.period(:,1)];
modper(:,2) = missing; modper(:,3) = S.period(:,2); modper(:,4)="start"; modper(:,5)=missing; modper(:,6)=S.period(:,3);
modcausd = [S.caus_decision(:,1)];
modcausd(:,2:3) = missing; modcausd(:,4)="caus_decision"; modcausd(:,5)=missing; modcausd(:,6)=S.caus_decision(:,3);
moddec = [cutDecision(:,1)];
moddec(:,2:3) = missing; moddec(:,4)="decision"; moddec(:,5)=cutDecision(:,2); moddec(:,6)=cutDecision(:,3);
modans = [Answers(:,1)];
modans(:,2:3) = missing; modans(:,4:6)=Answers(:,4:6); 
modcho = [S.choice(:,1)];
modcho(:,2:3) = missing; modcho(:,4)="choice"; modcho(:,5)=S.choice(:,2); modcho(:,6)=S.choice(:,3);

% trim end if col.2 is LOG
if modend(end,2) == "Log"
    modend(end,:) = [];
end

% pulling total logs
totalLogs = [modst; modend; modper; modcausd; moddec; modans; modcho]; 

%sort table
wholetable = splitvars(table(str2double(totalLogs(:,1:3)), totalLogs(:,4), str2double(totalLogs(:,5:6))));
wholetable.Properties.VariableNames = ["Time","Trial","Period","Event", "Note_direction", "Log"];
wholetable = sortrows(wholetable, [6 1 2 3]);
wholetable(:,2:3) = fillmissing(wholetable(:,2:3),'previous'); % fill missings


tablename = 'wholetable';

save(['C:\Users\sorin\Documents\MATLAB\23.03.16_Log error arrange\processed\' filenameo '\' filenameo '_' tablename], "wholetable");


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




























































