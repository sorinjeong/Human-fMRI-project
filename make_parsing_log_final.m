FileList = {'CL121121_1','CL121122_1','CL121128_1','CL121227_1','CL130107_1','CL130109_1','CL130114_2','CL130116_2',...
    'CL130121_2','CL130122_1','CL130130_1','CL130219_1','CL130220_1','CL130225_2','CL130226_1','CL130227_1'};


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
           S.(VarName{j})= [S.(VarName{j}); timeNevent.time(i), extractAfter(timeNevent.event(i), EventName{j})];   
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
save(['C:\Users\sorin\Documents\MATLAB\23.03.06_Log error arrange\processed\' filenameo], "S")

% writetable(timeNevent, "timeNevent.xlsx")

end