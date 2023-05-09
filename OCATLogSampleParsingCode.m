FileList = {'TaskDemo'};

%%
for fi = 1:numel(FileList)
    filenameo=FileList{fi};
    filefolder= ['D:\internship\MATLAB\23.05.09_O-CAT parsing\' filenameo];
    addpath(filefolder)

    %% Import the file
    Root =[ filefolder '\'];
    fileToRead1='2023.05.04-18.02.56.230Time_Behavior.csv';
    TimeBehavior = readcell(fileToRead1);

%% eventlogfile wo/coordinate
RawEventLog = [];
for i=1:height(TimeBehavior)
    if TimeBehavior{i,1}=="X"
        continue

    elseif TimeBehavior{i,1}=="LapStart"
        RawEventLog = [RawEventLog; NaN,TimeBehavior(i,1), TimeBehavior(i,4)];

    elseif TimeBehavior{i,3}=="Duration"
        RawEventLog = [RawEventLog; NaN, TimeBehavior(i,1),TimeBehavior(i,2); NaN, TimeBehavior(i,3), TimeBehavior(i,4)];

    elseif TimeBehavior{i,1}=="TrialType"
        RawEventLog = [RawEventLog; NaN,TimeBehavior(i,1), TimeBehavior(i,2)];

    else
        RawEventLog = [RawEventLog; TimeBehavior(i,4),TimeBehavior(i,1),TimeBehavior(i,3)];
    end
end

%% EventParsing

EventName = ["LapStart", "TrialStart", "TrialType", "ObjOn","ChoiceOn","Decision", "Duration", "ObjOff", "TrialEnd"];
VarName = ["Lap", "TrialStart", "Tr" + ...
    "ial", "Context", "Direction", "Location","Association", "Obj_ID", "ObjOn","ChoiceOn", "Choice","Decision", "Duration", "ObjOff", "TrialEnd"];

EventTable = table ;
Str=struct;
for i=1:height(RawEventLog)
    for k=1:length(VarName)
        if ~isfield(Str,VarName{k})
           Str.(VarName{k}) = [];
        end 

       % Lap
       if RawEventLog{i,2} == EventName(1)
           Str.Lap = [Str.Lap; RawEventLog{i,3}];
           break
       else


       % Trial
       if RawEventLog{i,2} == EventName(2)
           Str.TrialStart = [Str.TrialStart; RawEventLog{i,1}];
           if isempty(Str.Trial)
               Str.Trial = 1;
           else
               Str.Trial = [Str.Trial; Str.Trial(end)+1];
           end
           if Str.Trial(end) >= 5
               Str.Trial(end) = 1;
           end
           if isempty(Str.Lap) == 0
           Str.Lap = [Str.Lap; Str.Lap(end)];
           end
           break
       % TrialType
       elseif RawEventLog{i,2} == EventName(3)
           Numb = [];
           Numb = num2str(cell2mat(RawEventLog(i,3))) - '0';
           if Numb(1) == 1
               Str.Context = [Str.Context; 'F'];
           elseif Numb(1) == 2
               Str.Context = [Str.Context; 'C'];
           end
          Str.Direction = [Str.Direction; Numb(2)];
          Str.Location = [Str.Location; Numb(3)];
          Str.Association = [Str.Association; Numb(4)];
          Str.Obj_ID = [Str.Obj_ID; Numb(5)];
          if isempty(Str.Lap) == 0
           Str.Lap = [Str.Lap; Str.Lap(end)];
          end
          break
       %Object/choice On/Off, TrialEnd
       elseif RawEventLog{i,2} == EventName(4) | RawEventLog{i,2} == EventName(5) | RawEventLog{i,2} == EventName(8) | RawEventLog{i,2} == EventName(9)
           Str.(RawEventLog{i,2}) = [Str.(RawEventLog{i,2}); RawEventLog{i,1}];
          if isempty(Str.Lap) == 0
           Str.Lap = [Str.Lap; Str.Lap(end)];
          end
           break
       % decision, duration
       elseif RawEventLog{i,2} == EventName(6) | RawEventLog{i,2} == EventName(7)
           Str.(RawEventLog{i,2}) = [Str.(RawEventLog{i,2}); RawEventLog{i,3}];
           if isempty(Str.Lap) == 0
           Str.Lap = [Str.Lap; Str.Lap(end)];
          end
           break
       % Choice
       elseif contains(char(RawEventLog{i,2}),"Choice"+lettersPattern)
           Str.Choice = [Str.Choice; extractAfter(RawEventLog{i,2},"Choice")];
           if isempty(Str.Lap) == 0
           Str.Lap = [Str.Lap; Str.Lap(end)];
          end
           break
       end
           end
           break
    end
end


EventTable = table ;


          


end

