Subjects = {'Sub1', 'Sub2', 'Sub3', 'Sub4', 'Sub5'};

for fi = 1:numel(Subjects)
    subjectName=Subjects{fi};
    filefolder= ['D:\internship\MATLAB\23.06.02_O-CAT parsing_PVtask+TRlog\230602JiSunPilotLog'];
    cd(filefolder)

%% Import the file
    Root =[ filefolder '\'];
    fileToRead1= dir([subjectName '*Time_Behavior.csv']);
    TimeBehavior = readtable(fileToRead1.name);

%% Eventlogfile wo/coordinate 
    RawEventLog = TimeBehavior;
    Coord = find(contains(RawEventLog.Var1(:),'time'));
    RawEventLog(Coord,:) = [];
%% save TR log
    MREvent = find(contains(RawEventLog.Var1(:),'MR'));
    TRLog = RawEventLog(MREvent,:);

    if ~isfolder ([Root subjectName])
    mkdir([Root subjectName]); end
    save([Root subjectName '\' subjectName '_TRLog'], "TRLog");
    writetable(TRLog,[Root subjectName '\' subjectName '_TRLog.xlsx']);

    RawEventLog(MREvent,:) = [];

%% PV task Parsing
%pre
prePVTaskStart = find(contains(RawEventLog.Var1(:),'OCP_on'),1,"first");
prePVTaskEnd = find(contains(RawEventLog.Var1(:),'OCP_off'),1,"first");
prePVtaskLog = RawEventLog(prePVTaskStart:prePVTaskEnd,:);

save([Root subjectName '\' subjectName '_pre-PVtaskLog'], "prePVtaskLog");
writetable(prePVtaskLog,[Root subjectName '\' subjectName '_pre-PVtaskLog.xlsx']);
%post
postPVTaskStart = find(contains(RawEventLog.Var1(:),'OCP_on'),1,"last");
postPVTaskEnd = find(contains(RawEventLog.Var1(:),'OCP_off'),1,"last");
postPVtaskLog = RawEventLog(postPVTaskStart:postPVTaskEnd,:);

save([Root subjectName '\' subjectName '_post-PVtaskLog'], "postPVtaskLog");
writetable(postPVtaskLog,[Root subjectName '\' subjectName '_post-PVtaskLog.xlsx']);


%% MainTask Parsing
    VarName = ["Lap", "TrialStart", "Trial" "Context", "Direction", "Location","Association" ...
        , "Obj_ID", "ObjOn","ChoiceOn", "Choice","Decision", "Duration", "ObjOff", "TrialEnd"];

    %%Split Trial Type
    TrialTypeRow = find(contains(RawEventLog.Var3(:),'Type'));
    Numb = num2str(RawEventLog.Var4(TrialTypeRow)) - '0';

    ParsingPerTrial=struct;
    for s=1:5; ParsingPerTrial.(VarName{s+3}) = Numb(:,s); end
    for v=1:length(VarName); if ~isfield(ParsingPerTrial,VarName{v}); ParsingPerTrial.(VarName{v}) = [];end;end %event name별로 field 생성  
%% Event parsing
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

                ParsingPerTrial.Lap = [ParsingPerTrial.Lap; LapNumNTime(find(ParsingPerTrial.TrialStart(end) > LapNumNTime(:,1),1,'last'),2)];

            else
                % Choice
                if contains(EventName(i),"Choice"+("A"|"B"))
                    ParsingPerTrial.Choice = [ParsingPerTrial.Choice; string(extractAfter(EventName(i),"Choice"))];

                    % decision, duration
                elseif EventName(i)=="Decision"
                    ParsingPerTrial.Decision = [ParsingPerTrial.Decision; RawEventLog{i,2}];
                    ParsingPerTrial.Duration = [ParsingPerTrial.Duration; RawEventLog{i,4}];

                    %Object|choice On|Off, TrialEnd
                elseif ismember(EventName(i),VarName) & ~ismember(EventName(i),"Trial")
                    ParsingPerTrial.(EventName(i)) = [ParsingPerTrial.(EventName(i)); RawEventLog{i,2}];

                end
            end
            %choice missing (choiceOn과 objOff 이벤트 개수와 choice 개수가 다를 경우 missing 삽입)
            if length(ParsingPerTrial.ChoiceOn) == length(ParsingPerTrial.ObjOff) & length(ParsingPerTrial.ChoiceOn) ~= length(ParsingPerTrial.Choice)
                ParsingPerTrial.Choice = [ParsingPerTrial.Choice; missing];end           
    end

    % 1,2 --> forest(F), city (C)로 변경
    [r1,c1] = find(ParsingPerTrial.Context ==1);[r2,c2] = find(ParsingPerTrial.Context ==2); ParsingPerTrial.Context=string(ParsingPerTrial.Context);
    ParsingPerTrial.Context(r1,1)="F"; ParsingPerTrial.Context(r2,1)="C";



%% Save as Table 

% field 순서 맞추고 table화
ParsingPerTrial= orderfields(ParsingPerTrial,VarName);
LogTable=struct2table(ParsingPerTrial);

% save log data
save([Root subjectName '\' subjectName '_LogTable'], "LogTable");
writetable(LogTable,[Root subjectName '\' subjectName '_LogTable.xlsx']);


%% Plot
figure{fi}
clf
plottingX=1:height(LogTable);
hold on
yyaxis left
colororder({'#0072BD','#000000'})
title('RT & Correctness')
xlabel('Trial')
ylabel('RT')
RTplot = plot(LogTable,"Duration",'LineWidth',1.5,'Marker','.','MarkerSize',20);
Threshold = yline(1.5,'-.','Timeout','LabelHorizontalAlignment', 'center' ,'Color',"#0072BD");

xlim([1 height(LogTable)]);
ylim([0 2]);
pbaspect([2 1 1]);


yyaxis right

plottingY=LogTable.Decision';
[to_r, to_c] = find(plottingY == 2);
plottingY(1,to_c) = 1;


stairs(plottingY);
coloringX = [plottingX;plottingX];
coloringY = [plottingY;plottingY];
CorPlot = area(coloringX([2:end end]),coloringY(1:end));


TOmat = zeros(1,height(LogTable));
TOmat(1,to_c) = 1;
TOmat = [TOmat;TOmat];
TOPlot = area(coloringX([2:end end]),TOmat(1:end));
TOPlot.FaceColor = "#EDB120";
ylim([0 10]);
ylabel('correctness')
yticks([0 0.5 1])
yticklabels({'x', 'timeout', 'o'})
hold off

saveas(gcf,[Root subjectName '\' subjectName 'RTandCorrect.png'])

end


























