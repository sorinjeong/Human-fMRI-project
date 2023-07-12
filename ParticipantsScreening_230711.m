%% root 수기지정

cd('Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230711\');
load('Allsub_NumLogTable.mat');T=total_NumLogTable;clear("total_NumLogTable");
%Session info
SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
SubInfoFile = renamevars(SubInfoFile(2:19,[1 3 4 5]),["Var1","Var3","Var4","Var5"],["Session","PASS","Sex","Age"]);

%% session별 나누기
Subs = unique(T.Session,"rows","sorted");

SpliT=struct; To=table2array(T); To(To(:,12)==1,2:end)=missing;To=array2table(To);To.Properties.VariableNames = T.Properties.VariableNames;
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    SpliT.(varName) = To((To.Session(:) == Subs(i)),:);
% end; clear("varName", "i");


%% Context-Object 조합

C = SpliT.(varName).Context(find(SpliT.(varName).Association == 1));
O = SpliT.(varName).Object(find(SpliT.(varName).Association == 1));
A=[C O];
Asso= sort([A(find(O==4,1)),4; A(find(O==5,1)),5; A(find(O==6,1)),6; A(find(O==7,1)),7]);

N=num2str(Asso(1:4,2)');SubInfoFile.Combi_FFCC=zeros(height(SubInfoFile),1);
N(N ==' ')=[]; Combination = str2num(N);
SubInfoFile.Combi_FFCC(find(SubInfoFile.Session == Subs(i)))=Combination;

end; clear("varName", "i");

%% Data Group (Pass/Fail , Correct/Overall)

DataGroup = struct("PASS",[],"FAIL",[],"Correct",[],"Overall",[]);
DataGroup.Overall = SpliT;
P=SubInfoFile.Session(find(SubInfoFile.PASS == 1));
F=SubInfoFile.Session(find(SubInfoFile.PASS == 0));
for n=1:length(P); DataGroup.PASS.(sprintf('SUB_%.15g', P(n))) = SpliT.(sprintf('SUB_%.15g', P(n))); end
for n=1:length(F); DataGroup.FAIL.(sprintf('SUB_%.15g', F(n))) = SpliT.(sprintf('SUB_%.15g', F(n))); end

Tc=table2array(To); Tc(Tc(:,10)==0,2:end)=missing;Tc=array2table(Tc);Tc.Properties.VariableNames = To.Properties.VariableNames;
corr_Percent = [];halfbias=[];allbias=[];overall_corr = [];overall_sub =[];
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    DataGroup.Correct.(varName) = Tc((Tc.Session(:) == Subs(i)),:);

%% Performance-Bias
Screening.(varName)= struct("Accuracy",[],"RT",[],"Bias_all",[],"Bias_Lap",[],"Bias_Half",[]);
%Bias_overall trials
if find(Subs(i) == P)
ButtonA = length(find(DataGroup.PASS.(varName).Choice==1));
ButtonB = length(find(DataGroup.PASS.(varName).Choice==2));
Screening.(varName).Bias_all = (ButtonA-ButtonB) / height(DataGroup.PASS.(varName));
elseif find(Subs(i) == F)
ButtonA = length(find(DataGroup.FAIL.(varName).Choice==1));
ButtonB = length(find(DataGroup.FAIL.(varName).Choice==2));
Screening.(varName).Bias_all = (ButtonA-ButtonB) / height(DataGroup.FAIL.(varName));
end

%Bias_per lap
Screening.(varName).Bias_Lap = struct;
for j=1:8
    if find(Subs(i) == P)
        lapnum = find(DataGroup.PASS.(varName).Lap==j);
        lapnumchoice = DataGroup.PASS.(varName).Choice(lapnum);
    elseif find(Subs(i) == F)
        lapnum = find(DataGroup.FAIL.(varName).Lap==j);
        lapnumchoice = DataGroup.FAIL.(varName).Choice(lapnum);
    end

ButtonA_Lap = length(find(lapnumchoice==1));
ButtonB_Lap = length(find(lapnumchoice==2));
Screening.(varName).Bias_Lap.(sprintf('Lap%.15g', j)) = (ButtonA_Lap-ButtonB_Lap) / height(lapnumchoice);

%Bias_Last Half
if j==5
    if find(Subs(i) == P)
    halfchoice = DataGroup.PASS.(varName).Choice(lapnum(1):end);
    elseif find(Subs(i) == F)
    halfchoice = DataGroup.FAIL.(varName).Choice(lapnum(1):end);
    end

ButtonA_Half = length(find(halfchoice==1));
ButtonB_Half = length(find(halfchoice==2));
Screening.(varName).Bias_Half = (ButtonA_Half-ButtonB_Half) / height(halfchoice);
end
end

%% Performance-Accuracy
corr_Percent = [corr_Percent; ((height(DataGroup.Correct.(varName)))/32)*100];
allbias = [allbias; Screening.(varName).Bias_all];
halfbias = [halfbias; Screening.(varName).Bias_Half];


%% make correctness array
overall_sub = [overall_sub DataGroup.Overall.(varName).Session];

end
SubInfoFile = addvars(SubInfoFile, corr_Percent, allbias, halfbias); 



%% Screening_RT plot
% RT Y value 생성
overall_RT=[];D=struct;
for i=1:length(Subs)
    temp=T((T.Session(:) == Subs(i)),:);
    temp.RT(find(temp.isTimeout)) = nan;
    overall_RT = [overall_RT temp.RT] ;
end

% boxplot 생성
f=figure;
hold on
title('Response Time /Subject')
xlabel('Subject Number')
ylabel('RT(s)')

h = boxplot(overall_RT,SubInfoFile.Session);
idx = find(ismember([SubInfoFile.Session], F));
set(h(:,idx),'Color','red');
set(h(6,:),'Color','k');
% set(h(7,:),'MarkerEdgeColor','k');
set(h(1:2,:),'LineStyle','-');

ax = gca;
xTick = ax.XTick;
xLim = ax.XLim;
yLim = ax.YLim;

for i = 1:length(idx)
    text(xTick(idx(i)), yLim(1)-0.05*diff(yLim), ax.XTickLabel{idx(i)},...
        'Color', 'red', 'HorizontalAlignment', 'center','Rotation', 45);
ax.XTickLabel{idx(i)} = '';
end
meanValues = mean(overall_RT);
for i = 1:length(xTick)
    text(xTick(i), yLim(1)-0.1*diff(yLim), sprintf('Mean: %.2f',...
        meanValues(i)), 'HorizontalAlignment', 'center');
end































