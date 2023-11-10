%% 아래 두가지 variable 반드시 기입할 것!!!!

NumOfSubs = 1;
ParsingVersion = num2str(231109.2);

% %% subject numbering , folder root
% for subname = 1:NumOfSubs
%     Subjects{subname} = strcat("sub-",num2str(subname, '%02.f'));
% end
Subjects={"sub-01",	"sub-02",	"sub-03",	"sub-04",	"sub-05",	"sub-06",	"sub-07",	"sub-08",	"sub-09",	"sub-10"	,"sub-11",	"sub-12",	"sub-13",	"sub-14",	"sub-15",	"sub-16",	"sub-17",	"sub-18",	"sub-19",	"sub-20",	"sub-21.1","sub-21.2",	"sub-22",	"sub-23",	"sub-24",	"sub-25",	"sub-26",	"sub-27"};
% Subjects={"sub-21.1","sub-21.2",	"sub-22",	"sub-23",	"sub-24",	"sub-25",	"sub-26",	"sub-27"};
% Subjects={"sub-21.1"}

total_NumLogTable=[];
for fi = 1:numel(Subjects)
    Session=char(Subjects{fi});
    Root = ['Z:\E-Phys Analysis\fMRI_ocat\'];
    filefolder= [Root 'BehavLog\'];
    addpath(filefolder)
    cd(filefolder)

%add version
    savefolder= [Root 'BehavData_analyzed\ver_' ParsingVersion '\'];
    if ~isfolder (savefolder);mkdir(savefolder);mkdir([savefolder Session]);mkdir([savefolder 'GLM']);mkdir([savefolder 'Total']);end;addpath(savefolder)

%% Import the file
    fileToRead1= dir(fullfile(filefolder, '*.csv'));
    TimeBehavior = readtable(fileToRead1(fi).name);
    disp(['The current Session is: ' Session ', and the TimeBehavior table is: ' fileToRead1(fi).name]);  % Display the name of the Session and the TimeBehavior table
    SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
    SubInfoFile = renamevars(SubInfoFile(1:NumOfSubs,[1 3 4 5]),["SubjectNumber","HeadMovement_Attention","Sex","Age"],["Session","PASS","Sex","Age"]);

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

    disp([savefolder Session])
    if ~isfolder ([savefolder Session])
    mkdir([savefolder Session]);mkdir([savefolder 'GLM']);mkdir([savefolder 'Total']);end
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
    for s=1:5; ParsingPerTrial.(VarName{s+5}) = Numb(:,s); end
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
              
        % Correctness, RT, isTimeoutC
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
          if length(ParsingPerTrial.ChoiceOn) == length(ParsingPerTrial.ObjOff) & length(ParsingPerTrial.ChoiceOn) ~= length(ParsingPerTrial.Choice_txt)
                ParsingPerTrial.Choice_txt = [ParsingPerTrial.Choice_txt; missing];end   
       
 end

% Choice txt to Number
ParsingPerTrial.Choice_Num = NaN(height(ParsingPerTrial.Choice_txt), 1);  % Initialize Choice_Num with NaN
ParsingPerTrial.Choice_Num(ParsingPerTrial.Choice_txt == "A") = 1;
ParsingPerTrial.Choice_Num(ParsingPerTrial.Choice_txt == "B") = 2;

    %Context Num to txt
    ParsingPerTrial.Context_txt = string(ParsingPerTrial.Context_Num);
    ParsingPerTrial.Context_txt(find(ParsingPerTrial.Context_txt=="1"))="F";ParsingPerTrial.Context_txt(find(ParsingPerTrial.Context_txt=="2"))="C";

%% Make a Table
% field 순서 맞추고 table화
ParsingPerTrial= orderfields(ParsingPerTrial,VarName);
LogTable=struct2table(ParsingPerTrial);

%% Save the Table 

% save log data
save([savefolder Session '\' Session '_LogTable'], "LogTable");
writetable(LogTable,[savefolder Session '\' Session '_LogTable.xlsx']);

%save Total_Table
cd([savefolder 'Total']);
writetable(LogTable,'LogTable_Total.xlsx','Sheet', Session);
writetable(postPVtaskLog,'post-PV_Total.xlsx','Sheet', Session);
writetable(prePVtaskLog,'pre-PV_Total.xlsx','Sheet', Session);
writetable(TRLog,'TR_Total.xlsx','Sheet', Session);


%% Plot
cd([Root 'PilotData_analyzed\']) %legend 있는 곳으로 이동
clf
f=figure; f.Position; f.Position = [1500 500 1000 600];
%RT plot
plottingX=1:height(LogTable);
hold on
colororder({'#0072BD','#000000'})
yyaxis left
title([Session ': RT & Correctness'],"FontSize",18,"FontWeight","bold")


xlabel('Trial','FontSize',15,'FontWeight','bold')
ylabel('RT','FontSize',15,'FontWeight','bold')
RTplot = plot(LogTable,"RT",'Color','#0072BD', 'LineWidth',1.8,'Marker','.','MarkerSize',20);
Threshold = yline(1.5,'-.','Timeout','LabelHorizontalAlignment', 'center' ,'Color',"#0072BD",'LineWidth',1.2);
yticks(0:0.2:1.8)
xlim([1 height(LogTable)]);
ylim([-0.2 2]);
pbaspect([2 1 1]);

% correctness plot
yyaxis right

plottingY=LogTable.Correct_Num';
[to_r, to_c] = find(plottingY == 2);
plottingY(1,to_c) = 1;


stairs(plottingY);
coloringX = [plottingX;plottingX];
coloringY = [plottingY;plottingY];
CorPlot = area(coloringX([2:end end]),coloringY(1:end));
CorPlot.FaceColor = "k";
ylim([0 10]);

%legend
img = imread('Plot_legend.jpg');
image(img,'XData',[32 34],'YData',[3 0],'Clipping','off')

TOmat = zeros(1,height(LogTable));
TOmat(1,to_c) = 1;
TOmat = [TOmat;TOmat];
TOPlot = area(coloringX([2:end end]),TOmat(1:end));
TOPlot.FaceColor = "#EDB120";
yticks([])

hold off
saveas(gcf,[savefolder Session '\' Session '_Performance_Graph.png'])
close(f)

%% GLM corr/incorr
corr=[];incorr=[];
for co = 1:height(LogTable)
if ParsingPerTrial.Correct_Num(co) == 1
    corr= [corr ParsingPerTrial.ObjOn(co)];
else 
    incorr= [incorr ParsingPerTrial.ObjOn(co)];
end
end

save([savefolder 'GLM\corr_' Session] ,"corr",'-mat')
save([savefolder 'GLM\incorr_' Session] ,"incorr",'-mat')


%% Table for Analysis == LogTable_NumOnly

% "sub-" 문자열을 제거하고 숫자로 변환합니다.
subnum = str2double(strrep(Session, 'sub-', ''));
subnum = repmat(subnum, 32, 1);

LogTable_NumOnly=struct('Session',subnum,'Lap',LogTable.Lap,'Trial',LogTable.Trial,'Context',LogTable.Context_Num,'Direction',LogTable.Direction,...
    'Location',LogTable.Location,'Association',LogTable.Association,'Object',LogTable.Obj_ID,'Choice',LogTable.Choice_Num,...
    'Correct',LogTable.Correct_Num,'RT',LogTable.RT,'isTimeout',LogTable.isTimeout);



%% save / for 1 subject
LogTable_NumOnly=struct2table(LogTable_NumOnly);
save([savefolder Session '\' Session '_NumLogTable'], "LogTable_NumOnly");
writetable(LogTable_NumOnly,[savefolder Session '\' Session '_NumLogTable.xlsx']);

%% save /all subjects
total_NumLogTable=[total_NumLogTable; LogTable_NumOnly];
end
save([savefolder 'Allsub_NumLogTable'],"total_NumLogTable");
writetable(total_NumLogTable,[savefolder 'Allsub_NumLogTable.xlsx']);





% 
% %% Figures 모아놓기
% % 원본 디렉토리와 대상 디렉토리를 지정합니다.
% sourceDir = 'Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230828';
% targetDir = 'Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230828\Figures';  % 여기에 원하는 대상 디렉토리 경로를 입력하세요.
% 
% % 원본 디렉토리에서 모든 PNG 파일을 찾습니다.
% pngFiles = dir(fullfile(sourceDir, '**/*.png'));  % '**'는 하위 디렉토리를 포함하여 검색합니다.
% 
% % 각 PNG 파일을 대상 디렉토리로 복사합니다.
% for k = 1:length(pngFiles)
%     sourceFile = fullfile(pngFiles(k).folder, pngFiles(k).name);
%     copyfile(sourceFile, targetDir);
% end
% 
% disp('All PNG files have been copied to the target directory.');
% 
% 
% 
% %% figure 한번에 모은 figure 생성
% % 원본 디렉토리를 지정합니다.
% sourceDir = 'Z:\E-Phys Analysis\fMRI_ocat\BehavData_analyzed\ver_231109';
% 
% % 원본 디렉토리에서 모든 PNG 파일을 찾습니다.
% pngFiles = dir(fullfile(sourceDir, '**/*.png'));  % '**'는 하위 디렉토리를 포함하여 검색합니다.
% 
% % 새 figure를 생성합니다.
% figure;
% 
% % 각 PNG 파일을 subplot으로 표시합니다.
% numFiles = length(pngFiles);
% numRows = ceil(sqrt(numFiles));  % subplot의 행 수
% numCols = ceil(numFiles / numRows);  % subplot의 열 수
% for k = 1:numFiles
%     subplot(numRows, numCols, k);
%     img = imread(fullfile(pngFiles(k).folder, pngFiles(k).name));
%     imshow(img);
%     title(pngFiles(k).name);
% end
% 
% 


