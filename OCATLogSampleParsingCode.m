FileList = {'TaskDemo'};

for fi = 1:numel(FileList)
    filenameo=FileList{fi};
    filefolder= ['D:\internship\MATLAB\23.05.09_O-CAT parsing\' filenameo];
    addpath(filefolder)

Import the file
    Root =[ filefolder '\'];
    fileToRead1='2023.05.04-18.02.56.230Time_Behavior.csv';
    TimeBehavior = readtable(fileToRead1);

Eventlogfile wo/coordinate
    RawEventLog = TimeBehavior;
    Coord = find(contains(RawEventLog.Var1(:),'X'));
    RawEventLog(Coord,:) = [];

    VarName = ["Lap", "TrialStart", "Trial", "Context", "Direction", "Location","Association" ...
        , "Obj_ID", "ObjOn","ChoiceOn", "Choice","Decision", "Duration", "ObjOff", "TrialEnd"];

    %% Split Trial Type
    TrialTypeRow = find(contains(RawEventLog.Var1(:),'TrialType'));
    Numb = num2str(RawEventLog.Var2(TrialTypeRow)) - '0';

    ParsingPerTrial=struct;
    for s=1:5; ParsingPerTrial.(VarName{s+3}) = Numb(:,s); end
    for v=1:length(VarName); if ~isfield(ParsingPerTrial,VarName{v}); ParsingPerTrial.(VarName{v}) = [];end;end %event name별로 field 생성
 Event parsing
    LapNumNTime = [];EventName=[];EventName = string(RawEventLog.Var1(:));
    for i=1:height(RawEventLog)
        if EventName(i)=="LapStart"
            LapNumNTime = [LapNumNTime; RawEventLog.Var2(i), RawEventLog.Var4(i)];
        end
              % Lap, Trial
            if EventName(i)=="TrialStart"
                ParsingPerTrial.TrialStart = [ParsingPerTrial.TrialStart; RawEventLog{i,2}];

                if isempty(ParsingPerTrial.Trial) | ParsingPerTrial.Trial(end)==4; t=1;
                elseif ParsingPerTrial.Trial(end)==1; t=2;elseif ParsingPerTrial.Trial(end)==2; t=3;elseif ParsingPerTrial.Trial(end)==3; t=4;end
                ParsingPerTrial.Trial = [ParsingPerTrial.Trial; t];

                ParsingPerTrial.Lap = [ParsingPerTrial.Lap; LapNumNTime(find(ParsingPerTrial.TrialStart(end) > LapNumNTime(:,1),1,'last'),2)]

            else
                % Choice
                if contains(EventName(i),"Choice"+("A"|"B"))
                    ParsingPerTrial.Choice = [ParsingPerTrial.Choice; string(extractAfter(EventName(i),"Choice"))];

                    % decision, duration
                elseif EventName(i)=="Decision"
                    ParsingPerTrial.Decision = [ParsingPerTrial.Decision; RawEventLog{i,2}];
                    ParsingPerTrial.Duration = [ParsingPerTrial.Duration; RawEventLog{i,4}];

                    %Object|choice On|Off, TrialEnd
                elseif ismember(EventName(i),VarName)
                    ParsingPerTrial.(EventName(i)) = [ParsingPerTrial.(EventName(i)); RawEventLog{i,2}];

                end
            end
            %choice missing (choiceOn과 objOff 이벤트 개수와 choice 개수가 다를 경우 missing 삽입)
            if length(ParsingPerTrial.ChoiceOn) == length(ParsingPerTrial.ObjOff) & length(ParsingPerTrial.ChoiceOn) ~= length(ParsingPerTrial.Choice)
                ParsingPerTrial.Choice = [ParsingPerTrial.Choice; missing];end           
    end

    % 1,2 --> forest(F), city (C)로 변경
    [r1,c1] = find(ParsingPerTrial.Context ==1);[r2,c2] = find(ParsingPerTrial.Context ==2); ParsingPerTrial.Context=string(ParsingPerTrial.Context)
    ParsingPerTrial.Context(r1,1)="F"; ParsingPerTrial.Context(r2,1)="C"

end

% field 순서 맞추고 table화
ParsingPerTrial= orderfields(ParsingPerTrial,VarName)
LogTable=struct2table(ParsingPerTrial)


%% save log data

if ~isfolder (['D:\internship\MATLAB\23.05.09_O-CAT parsing\processed\' filenameo])
mkdir(['D:\internship\MATLAB\23.05.09_O-CAT parsing\processed\' filenameo]); end

save(['D:\internship\MATLAB\23.05.09_O-CAT parsing\processed\' filenameo '\' filenameo], "LogTable");


