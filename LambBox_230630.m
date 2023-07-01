%% subject numbering , folder root
for subname = 1:1
    Subjects{subname} = sprintf('Sub%.15g', subname);
end
% for fi = 1:numel(Subjects)
%     subjectName=Subjects{fi};
%     Root = ['Z:\E-Phys Analysis\fMRI_ocat\'];
%     filefolder= [Root 'PilotData'];
%     addpath(filefolder)
%     cd(filefolder)
%     % savefolder= [Root 'PilotData_analyzed\'];
%     savefolder='D:\internship\MATLAB\23.06.30_LambBox\';
%     addpath(savefolder)

for fi = 1:numel(Subjects)
    subjectName=Subjects{fi};
    % Root = ['Z:\E-Phys Analysis\fMRI_ocat\'];
   Root = ['D:\internship\MATLAB\23.06.30_LambBox'];
    % filefolder= [Root 'PilotData'];
    filefolder=Root;
    addpath(filefolder)
    cd(filefolder)
    % savefolder= [Root 'PilotData_analyzed\'];
    savefolder= [Root '\'];
    addpath(savefolder)


%% Import the file
    fileToRead1= dir([subjectName '*Time_Behavior.csv']);
    TimeBehavior = readtable(fileToRead1.name);

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

    ParsingPerTrial=struct;
    for s=1:5; ParsingPerTrial.(VarName{s+3}) = Numb(:,s); end
    for v=1:length(VarName); if ~isfield(ParsingPerTrial,VarName{v}); ParsingPerTrial.(VarName{v}) = [];end;end %event name별로 field 생성  
%% Event parsing
    LapNumNTime = [];EventName=[];EventName = string(RawEventLog.Var1(:));Trial1to32=[];
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
                    if RawEventLog{i,4} > 1.5 & RawEventLog{i,2} ~= 2
                        RawEventLog{i,2} = 2;end
                    ParsingPerTrial.Decision = [ParsingPerTrial.Decision; RawEventLog{i,2}];
                    ParsingPerTrial.Duration = [ParsingPerTrial.Duration; RawEventLog{i,4}];

                    %Object|choice On|Off, TrialEnd
                elseif ismember(EventName(i),VarName) & ~ismember(EventName(i),"Trial")
                    ParsingPerTrial.(EventName(i)) = [ParsingPerTrial.(EventName(i)); RawEventLog{i,2}];

                end
                %trial 1 to 32
            if EventName(i)=="Trial"
                Trial1to32 = [Trial1to32; RawEventLog{i,2}];end
            end
            %choice missing (choiceOn과 objOff 이벤트 개수와 choice 개수가 다를 경우 missing 삽입)
            if length(ParsingPerTrial.ChoiceOn) == length(ParsingPerTrial.ObjOff) & length(ParsingPerTrial.ChoiceOn) ~= length(ParsingPerTrial.Choice)
                ParsingPerTrial.Choice = [ParsingPerTrial.Choice; missing];end           
    end

    % 1,2 --> forest(F), city (C)로 변경
    ContextNum=ParsingPerTrial.Context;
    [r1,c1] = find(ParsingPerTrial.Context ==1);[r2,c2] = find(ParsingPerTrial.Context ==2); ParsingPerTrial.Context=string(ParsingPerTrial.Context);
    ParsingPerTrial.Context(r1,1)="F"; ParsingPerTrial.Context(r2,1)="C";


%% Make a Table
% field 순서 맞추고 table화
ParsingPerTrial= orderfields(ParsingPerTrial,VarName);
LogTable=struct2table(ParsingPerTrial);

% %Make Bias Table
% for l=1:height(LogTable)
%    for lap=1:4:32
%         if l==lap & l~=1
%             LeftButton=0;RightButton=0;
%         if contains(LogTable.Choice(l),"A"); LeftButton = LeftButton+1; else; RightButton = RightButton+1;end
%         else
%             if contains(LogTable.Choice(l),"A"); LeftButton = LeftButton+1; else; RightButton = RightButton+1;end
%         LogTable.Bias(lap) = (LeftButton - RightButton)/4;end;end;end
%         % end;
%        % if contains(LogTable.Choice(l),"A"); LeftButton = LeftButton+1; else; RightButton = RightButton+1;end;end;end


%% Save the Table 
% save log data
save([savefolder subjectName '\' subjectName '_LogTable'], "LogTable");
writetable(LogTable,[savefolder subjectName '\' subjectName '_LogTable.xlsx']);
%save trial 1 to 32, context number
save([savefolder subjectName '\' subjectName '_Trial1to32'], "Trial1to32");
save([savefolder subjectName '\' subjectName '_ContextNum'], "ContextNum");

%save Total_Table
cd([Root '\analyzed']);
writetable(LogTable,'TotalSubject_LogTable.xlsx','Sheet', subjectName);
writetable(postPVtaskLog,'TotalSubject_post-PV.xlsx','Sheet', subjectName);
writetable(prePVtaskLog,'TotalSubject_pre-PV.xlsx','Sheet', subjectName);
writetable(TRLog,'TotalSubject_TR.xlsx','Sheet', subjectName);

%% Plot
cd(savefolder) %legend 있는 곳으로 이동
clf
f=figure; f.Position; f.Position = [1500 500 1000 600];
%RT plot
plottingX=1:height(LogTable);
hold on
colororder({'#0072BD','#000000'})
yyaxis left
title([subjectName ': RT & Correctness'],"FontSize",15,"FontWeight","bold")


xlabel('Trial','FontSize',13,'FontWeight','bold')
ylabel('RT','FontSize',13,'FontWeight','bold')
RTplot = plot(LogTable,"Duration",'Color','#0072BD', 'LineWidth',1.5,'Marker','.','MarkerSize',20);
Threshold = yline(1.5,'-.','Timeout','LabelHorizontalAlignment', 'center' ,'Color',"#0072BD",'LineWidth',1.2);
yticks(0:0.2:1.8)
xlim([1 height(LogTable)]);
ylim([-0.2 2]);
pbaspect([2 1 1]);

% correctness plot
yyaxis right

plottingY=LogTable.Decision';
[to_r, to_c] = find(plottingY == 2);
plottingY(1,to_c) = 1;


stairs(plottingY);
coloringX = [plottingX;plottingX];
coloringY = [plottingY;plottingY];
CorPlot = area(coloringX([2:end end]),coloringY(1:end));
CorPlot.FaceColor = "k";
ylim([0 10]);

%legend
img = imread('Plot_legend_v3.jpg');
image(img,'XData',[32 36],'YData',[3 -0.5],'Clipping','off')

TOmat = zeros(1,height(LogTable));
TOmat(1,to_c) = 1;
TOmat = [TOmat;TOmat];
TOPlot = area(coloringX([2:end end]),TOmat(1:end));
TOPlot.FaceColor = "#EDB120";
yticks([])

hold off

saveas(gcf,[savefolder subjectName '\' subjectName '_Performance_Graph.png'])
close(f)

%% GLM corr/incorr
corr=[];incorr=[];
for co = 1:height(LogTable)
if ParsingPerTrial.Decision(co) == 1
    corr= [corr ParsingPerTrial.ObjOn(co)];
else 
    incorr= [incorr ParsingPerTrial.ObjOn(co)];
end
end

save([savefolder 'GLM\corr_' subjectName] ,"corr",'-mat')
save([savefolder 'GLM\incorr_' subjectName] ,"incorr",'-mat')
end




%% LambBox for Analysis

%column: session, ocntext, location, object, assoc, choice, corr, RT,
%Timeout, trial, lap





















