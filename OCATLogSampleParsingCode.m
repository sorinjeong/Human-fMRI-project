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

%% EventNameVariable
EventName = ["LapStart", "TrialStart", "TrialType", "ObjOn","ChoiceOn","Decision", "Duration", "ObjOff", "TrialEnd"];

%% struct wo/coordinate
EventLog = [];
for i=1:height(TimeBehavior)
    if TimeBehavior{i,1}=="X"
        continue

    elseif TimeBehavior{i,1}=="LapStart"
        EventLog = [EventLog; NaN,TimeBehavior(i,1), TimeBehavior(i,4)];

    elseif TimeBehavior{i,3}=="Duration"
        EventLog = [EventLog; NaN, TimeBehavior(i,1),TimeBehavior(i,2); NaN, TimeBehavior(i,3), TimeBehavior(i,4)];

    elseif TimeBehavior{i,1}=="TrialType"
        EventLog = [EventLog; NaN,TimeBehavior(i,1), TimeBehavior(i,2)];

    else
        EventLog = [EventLog; TimeBehavior(i,4),TimeBehavior(i,1),TimeBehavior(i,3)];
    end
end





































end

