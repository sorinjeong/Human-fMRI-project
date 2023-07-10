NumOfSubs = 18;

%% subject numbering , folder root
for subname = 1:NumOfSubs
    Subjects{subname} = sprintf('Sub%.15g', (subname + 84));
end
total_lambBox=[];
for fi = 1:numel(Subjects)
    Session=Subjects{fi};
    Root = ['Z:\E-Phys Analysis\fMRI_ocat\'];
    filefolder= [Root 'PilotData\'];
    addpath(filefolder)
    cd(filefolder)
    savefolder= [Root 'PilotData_analyzed\'];
    addpath(savefolder)


%% Import the file
    fileToRead1= dir([Session '_*Time_Behavior.csv']);
    TimeBehavior = readtable(fileToRead1.name);
    SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
    SubInfoFile = renamevars(SubInfoFile(2:NumOfSubs+1,[1 3 4]),["Var1","Var3","Var4"],["Session","Sex","Age"]);

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


    if ~isfolder ([savefolder Session])
    mkdir([savefolder Session]); end
    save([savefolder Session '\' Session '_TRLog'], "TRLog");
    writetable(TRLog,[savefolder Session '\' Session '_TRLog.xlsx']);

    RawEventLog(MREvent,:) = [];


%% PV task Parsing
%pre
prePVTaskStart = find(contains(RawEventLog.Var1(:),'OCP_on'),1,"first");
prePVTaskEnd = find(contains(RawEventLog.Var1(:),'OCP_off'),1,"first");
prePVtaskLog = RawEventLog(prePVTaskStart:prePVTaskEnd,:);

save([savefolder Session '\' Session '_pre-PVtaskLog'], "prePVtaskLog");
writetable(prePVtaskLog,[savefolder Session '\' Session '_pre-PVtaskLog.xlsx']);
%post
postPVTaskStart = find(contains(RawEventLog.Var1(:),'OCP_on'),1,"last");
postPVTaskEnd = find(contains(RawEventLog.Var1(:),'OCP_off'),1,"last");
postPVtaskLog = RawEventLog(postPVTaskStart:postPVTaskEnd,:);

save([savefolder Session '\' Session '_post-PVtaskLog'], "postPVtaskLog");
writetable(postPVtaskLog,[savefolder Session '\' Session '_post-PVtaskLog.xlsx']);


%% MainTask Parsing
    VarName = ["Lap", "TrialStart", "Trial","Lap_Trial", "Context_txt", "Context_Num", "Direction", "Location","Association" ...
        , "Obj_ID", "ObjOn","ChoiceOn", "Choice_Num","Choice_txt","Correct_Num","Correct_txt","isTimeout" "RT", "ObjOff", "TrialEnd"];

    %%Split Trial Type
    TrialTypeRow = find(contains(RawEventLog.Var3(:),'Type'));
    Numb = num2str(RawEventLog.Var4(TrialTypeRow)) - '0';
    %event name별로 field 생성  
    ParsingPerTrial=struct;
    for s=1:5; ParsingPerTrial.(VarName{s+4}) = Numb(:,s); end
    for v=1:length(VarName); if ~isfield(ParsingPerTrial,VarName{v}); ParsingPerTrial.(VarName{v}) = [];end;end 

    %% Event Parsing
 LapNumNTime = [];EventName=[];EventName = string(RawEventLog.Var1(:));
 for i=1:height(RawEventLog)
      % Lap, Trial
     if EventName(i)=="LapStart"
         LapNumNTime = [LapNumNTime; RawEventLog.Var2(i), RawEventLog.Var4(i)];end
     if EventName(i)=="TrialStart"
        ParsingPerTrial.TrialStart = [ParsingPerTrial.TrialStart; RawEventLog{i,2}];
     elseif EventName(i)=="Trial"
        ParsingPerTrial.Trial = [ParsingPerTrial.Trial; RawEventLog{i,2}];
            %Lap_Trial
            if isempty(ParsingPerTrial.Lap_Trial) | ParsingPerTrial.Lap_Trial(end)==4; t=1;
            elseif ParsingPerTrial.Lap_Trial(end)==1; t=2;elseif ParsingPerTrial.Lap_Trial(end)==2; t=3;elseif ParsingPerTrial.Lap_Trial(end)==3; t=4;end
        ParsingPerTrial.Lap_Trial = [ParsingPerTrial.Lap_Trial; t];
            %Lap
        ParsingPerTrial.Lap = [ParsingPerTrial.Lap; LapNumNTime(find(ParsingPerTrial.TrialStart(end) > LapNumNTime(:,1),1,'last'),2)];
       
     else
         % Choice
             if contains(EventName(i),"Choice"+("A"|"B"))
              ParsingPerTrial.Choice_txt = [ParsingPerTrial.Choice_txt; string(extractAfter(EventName(i),"Choice"))];
            
        % Correctness, RT, isTimeout
             elseif EventName(i)=="Decision"
                    if RawEventLog{i,4} > 1.5 & RawEventLog{i,2} ~= 2
                    RawEventLog{i,2} = 2;end
              ParsingPerTrial.Correct_Num = [ParsingPerTrial.Correct_Num; RawEventLog{i,2}];
              ParsingPerTrial.RT = [ParsingPerTrial.RT; RawEventLog{i,4}];
                %txt, timeout
                   if RawEventLog{i,2} == 0; ParsingPerTrial.Correct_txt = [ParsingPerTrial.Correct_txt; "Incorrect"];
                                          ParsingPerTrial.isTimeout = [ParsingPerTrial.isTimeout; 0];
                  elseif RawEventLog{i,2} == 1; ParsingPerTrial.Correct_txt = [ParsingPerTrial.Correct_txt; "Correct"];
                                              ParsingPerTrial.isTimeout = [ParsingPerTrial.isTimeout; 0];
                  elseif RawEventLog{i,2} == 2; ParsingPerTrial.Correct_txt = [ParsingPerTrial.Correct_txt; "TimeOut"];
                                              ParsingPerTrial.isTimeout = [ParsingPerTrial.isTimeout; 1];
                   end

                
              elseif ismember(EventName(i),VarName) & ~ismember(EventName(i),"Trial")
                  ParsingPerTrial.(EventName(i)) = [ParsingPerTrial.(EventName(i)); RawEventLog{i,2}];
              end
        end

        %choice missing (choiceOn과 objOff 이벤트 개수와 choice 개수가 다를 경우 missing 삽입)
          if length(ParsingPerTrial.ChoiceOn) == length(ParsingPerTrial.ObjOff) & length(ParsingPerTrial.ChoiceOn) ~= length(ParsingPerTrial.Choice)
                ParsingPerTrial.Choice_txt = [ParsingPerTrial.Choice_txt; missing];end   
        
         ParsingPerTrial.Choice_Num(find(ParsingPerTrial.Choice_txt=="A"))=1;ParsingPerTrial.Choice_Num(find(ParsingPerTrial.Choice_txt=="B"))=2;
         ParsingPerTrial.Choice_Num(find(ParsingPerTrial.Choice_txt==missing))=missing;ParsingPerTrial.Choice_Num=ParsingPerTrial.Choice_Num';
     end

     %Context Num to txt
    ParsingPerTrial.Context_txt(find(ParsingPerTrial.Context_Num==1))="F";ParsingPerTrial.Context_txt(find(ParsingPerTrial.Context_Num==2))="C";
    ParsingPerTrial.Context_txt=ParsingPerTrial.Context_txt';


%% Make a Table
% field 순서 맞추고 table화
ParsingPerTrial= orderfields(ParsingPerTrial,VarName);
LogTable=struct2table(ParsingPerTrial);

%% Save the Table 

% save log data
save([savefolder subjectName '\' subjectName '_LogTable'], "LogTable");
writetable(LogTable,[savefolder subjectName '\' subjectName '_LogTable.xlsx']);
%save trial 1 to 32, context number
save([savefolder subjectName '\' subjectName '_Trial1to32'], "Trial1to32");
save([savefolder subjectName '\' subjectName '_ContextNum'], "ContextNum");

%save Total_Table
cd(savefolder);
writetable(LogTable,'TotalSubject_LogTable.xlsx','Sheet', subjectName);
writetable(postPVtaskLog,'TotalSubject_post-PV.xlsx','Sheet', subjectName);
writetable(prePVtaskLog,'TotalSubject_pre-PV.xlsx','Sheet', subjectName);
writetable(TRLog,'TotalSubject_TR.xlsx','Sheet', subjectName);


































