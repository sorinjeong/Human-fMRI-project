NumOfSubs = 18;




%% subject numbering , folder root
for subname = 1:NumOfSubs
    Subjects{subname} = sprintf('Sub%.15g', (subname + 84));
end
total_lambBox=[];
for fi = 1:numel(Subjects)
    subjectName=Subjects{fi};
    Root = ['Z:\E-Phys Analysis\fMRI_ocat\'];
    filefolder= [Root 'PilotData\'];
    addpath(filefolder)
    cd(filefolder)
    savefolder= [Root 'PilotData_analyzed\'];
    addpath(savefolder)


%% Import the file
    fileToRead1= dir([subjectName '_*Time_Behavior.csv']);
    TimeBehavior = readtable(fileToRead1.name);
    SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
    SubInfoFile = SubInfoFile(1:length(Subjects),[1 3 4]);

%% Eventlogfile wo/coordinate 
    RawEventLog = TimeBehavior;
    Coord = find(contains(RawEventLog.Var1(:),'time'));
    RawEventLog(Coord,:) = [];
%% save TR log
    MREvent = find(contains(RawEventLog.Var1(:),'MR'));
    MRStartTime = RawEventLog.Var2(MREvent(1));
    % Time - MRStartTime
    for i=1:height(RawEventLog)
        if RawEventLog.Var1(i) ~= "Decision" & RawEventLog.Var1(i) ~= "Trial"
            RawEventLog.Var2(i) =  minus(RawEventLog.Var2(i), MRStartTime);
        else 
            continue
        end
    end
    TRLog = RawEventLog(MREvent,:);


    if ~isfolder ([savefolder subjectName])
    mkdir([savefolder subjectName]); end
    save([savefolder subjectName '\' subjectName '_TRLog'], "TRLog");
    writetable(TRLog,[savefolder subjectName '\' subjectName '_TRLog.xlsx']);

    RawEventLog(MREvent,:) = [];


%% PV task Parsing
%pre
prePVTaskStart = find(contains(RawEventLog.Var1(:),'OCP_on'),1,"first");
prePVTaskEnd = find(contains(RawEventLog.Var1(:),'OCP_off'),1,"first");
prePVtaskLog = RawEventLog(prePVTaskStart:prePVTaskEnd,:);

save([savefolder subjectName '\' subjectName '_pre-PVtaskLog'], "prePVtaskLog");
writetable(prePVtaskLog,[savefolder subjectName '\' subjectName '_pre-PVtaskLog.xlsx']);
%post
postPVTaskStart = find(contains(RawEventLog.Var1(:),'OCP_on'),1,"last");
postPVTaskEnd = find(contains(RawEventLog.Var1(:),'OCP_off'),1,"last");
postPVtaskLog = RawEventLog(postPVTaskStart:postPVTaskEnd,:);

save([savefolder subjectName '\' subjectName '_post-PVtaskLog'], "postPVtaskLog");
writetable(postPVtaskLog,[savefolder subjectName '\' subjectName '_post-PVtaskLog.xlsx']);


%% MainTask Parsing
    VarName = ["Lap", "TrialStart", "Trial" "Context", "Direction", "Location","Association" ...
        , "Obj_ID", "ObjOn","ChoiceOn", "Choice","Decision", "Duration", "ObjOff", "TrialEnd"];

    %%Split Trial Type
    TrialTypeRow = find(contains(RawEventLog.Var3(:),'Type'));
    Numb = num2str(RawEventLog.Var4(TrialTypeRow)) - '0';

%% Event Parsing
 LapNumNTime = [];EventName=[];EventName = string(RawEventLog.Var1(:));
 for i=1:height(RawEventLog)
     if EventName(i)=="LapStart"
         LapNumNTime = [LapNumNTime; RawEventLog.Var2(i), RawEventLog.Var4(i)];
        end
              % Lap, Trial
            if EventName(i)=="TrialStart"












































