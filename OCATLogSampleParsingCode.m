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
    LapNumNTime = [];


    for i=1:height(RawEventLog)
        if RawEventLog.Var1(i)=="LapStart"
            LapNumNTime = [LapNumNTime; RawEventLog.Var2(i), RawEventLog.Var4(i)];
        end

     
            if RawEventLog.Var1(i)=="TrialStart"
                ParsingPerTrial.TrialStart = [ParsingPerTrial.TrialStart; RawEventLog{i,2}];
                ParsingPerTrial.Trial = [ParsingPerTrial.Trial; RawEventLog{i,4}];

            else
                for v=1:length(VarName)
                    if RawEventLog.Var1(i) == VarName(v)
                        name = VarName(v)
                % Choice
                if contains(char(RawEventLog{i,1}),"Choice"+('A'|'B'))
                    ParsingPerTrial.(name) = [ParsingPerTrial.(name); extractAfter(char(RawEventLog.Var1(i)),"Choice")];

                    % decision, duration
                elseif RawEventLog.Var1(i)=="Decision"
                    ParsingPerTrial.Decision = [ParsingPerTrial.Decision; RawEventLog{i,2}];
                    ParsingPerTrial.Duration = [ParsingPerTrial.Decision; RawEventLog{i,4}];

                    %Object|choice On|Off, TrialEnd
                else
                    ParsingPerTrial.(name) = [ParsingPerTrial.(name); RawEventLog{i,2}];
                end
                    end
                end
            end
        end

    for j=1:height(LapNumNTime)
        if isempty(ParsingPerTrial.Lap)

            ParsingPerTrial.Lap = LapNumNTime(1,2);
        elseif ParsingPerTrial.Trial(end) > LapNumNTime(j,1)
            ParsingPerTrial.Lap = [ParsingPerTrial.Lap; LapNumNTime(j,2)];
        elseif ParsingPerTrial.Trial(end) < LapNumNTime(j,1)
            continue
        end
    end



end



