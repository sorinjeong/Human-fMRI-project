%% root 수기지정

cd('Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230711\');
load('Allsub_NumLogTable.mat');T=total_NumLogTable;clear("total_NumLogTable");
%Session info
SubInfoFile = readtable('OCAT subject info (pilot).xlsx','ReadRowNames',false);
SubInfoFile = renamevars(SubInfoFile(2:19,[1 3 4 5]),["Var1","Var3","Var4","Var5"],["Session","PASS","Sex","Age"]);

%% session별 나누기
Subs = unique(T.Session,"rows","sorted");

SpliT=struct; To=table2array(T); To(To(:,12)==1,2:end)=missing;To=array2table(To);To.Properties.VariableNames = T.Properties.VariableNames;
Combination=[];
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    SpliT.(varName) = To((To.Session(:) == Subs(i)),:);
% end; clear("varName", "i");


%% Context-Object 조합

C = SpliT.(varName).Context(find(SpliT.(varName).Association == 1));
O = SpliT.(varName).Object(find(SpliT.(varName).Association == 1));
A=[C O];
Asso= sortrows([A(find(O==4,1)),4; A(find(O==5,1)),5; A(find(O==6,1)),6; A(find(O==7,1)),7]);

N=num2str(Asso(1:4,2)');SubInfoFile.Combi_FFCC=zeros(height(SubInfoFile),1);
N(N ==' ')=[]; Combination = [Combination; str2double(N)];
end; SubInfoFile.Combi_FFCC=Combination; clear("varName", "i");

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
corr_Percent = [corr_Percent; ((length(find(DataGroup.Correct.(varName).Correct==1)))/32)*100];
allbias = [allbias; Screening.(varName).Bias_all];
halfbias = [halfbias; Screening.(varName).Bias_Half];

%% make correctness array
overall_sub = [overall_sub DataGroup.Overall.(varName).Session];

end
SubInfoFile = addvars(SubInfoFile, corr_Percent, allbias, halfbias); 



%% Screening_RT plot - for each subject
% overall_RT, correct trial에서의 RT 생성
overall_RT=[];corr_RT=[];incorr_RT=[];
for i=1:length(Subs)
    temp=T((T.Session(:) == Subs(i)),:);
    temp.RT(find(temp.isTimeout)) = nan;
    overall_RT = [overall_RT temp.RT] ;
    ttemp=temp;
    temp.RT(find(temp.Correct==0)) = nan;
    corr_RT = [corr_RT temp.RT] ;
    ttemp.RT(find(ttemp.Correct)) = nan;
    incorr_RT = [incorr_RT ttemp.RT];
end
clear("temp")

%screening RT mean값 구해서 넣고 accuracy도 채워넣기
meanValues = nanmean(overall_RT);fns = fieldnames(Screening);
for i=1:length(Subs)
    Screening.(fns{i}).RT = meanValues(i);
    Screening.(fns{i}).Accuracy = corr_Percent(i);
end

%% boxplot 생성
f=figure;f.Position = [1645,857,829,594]; f.PaperSize = [21,30];
hold on
title('Response Time (for each Subject)')
xlabel('Subject Number')
ylabel('RT(s)')

h = boxplot(overall_RT,SubInfoFile.Session, OutlierSize=10^(-200));
fail_group = find(ismember([SubInfoFile.Session], F));
set(h(:,fail_group),'Color','red');
set(h(6,:),'Color','k');
set(h(7,:),'MarkerEdgeColor','w');
set(h(1:2,:),'LineStyle','-');

ax = gca;
xTick = ax.XTick;
xLim = ax.XLim;
ylim([0 1.6]);
yLim = ylim;

for i = 1:length(fail_group)
    text(xTick(fail_group(i)), yLim(1)-0.03*diff(yLim), ax.XTickLabel{fail_group(i)},...
        'Color', 'red', 'HorizontalAlignment', 'center');
ax.XTickLabel{fail_group(i)} = '';
end

% % plot(meanValues) % mean값 넣고싶다면 이 주석을 풀어!
% meanValues = nanmean(overall_RT);
% for i = 1:18
%     text(i, 1.55, jjnum2str(meanValues(i),2), 'HorizontalAlignment', 'center');
% end


%% Screening_RT plot - Corr/inCorr

% boxplot 생성 - all sub, corr/incorr 2 boxes
f=figure;f.Position = [1645,857,829,594]; f.PaperSize = [21,30];
hold on
title('Response Time (for Correct and Incorrect Trials)')
xlabel('Correctness')
ylabel('RT(s)')

group = [ones(1,18), 2*ones(1,18)];
labels = {'Correct', 'Incorrect'};
xlim([0 3]);ylim([0 1.3]);
ha = boxplot([nanmean(corr_RT) nanmean(incorr_RT)],group, 'Labels', labels);


% boxplot 생성 - PASS sub, corr/incorr 2boxes

f=figure;f.Position = [1645,857,829,594]; f.PaperSize = [21,30];
hold on
title('Response Time (for Correct and Incorrect Trials) by Group(P/F)')
xlabel('Correctness')
ylabel('RT(s)')
xlim([0 5]);ylim([0 1.3]);
pass_group = find(ismember([SubInfoFile.Session], P));

% 각 그룹의 데이터 계산하기
pass_corr_RT = nanmean(corr_RT(:, pass_group));
pass_incorr_RT = nanmean(incorr_RT(:, pass_group));
fail_corr_RT = nanmean(corr_RT(:, fail_group));
fail_incorr_RT = nanmean(incorr_RT(:, fail_group));

% boxplot 그리기
group = [ones(1,length(pass_corr_RT)), 2*ones(1,length(pass_incorr_RT)), ...
         3*ones(1,length(fail_corr_RT)), 4*ones(1,length(fail_incorr_RT))];
labels = {'Correct (Pass)', 'Incorrect (Pass)', 'Correct (Fail)', 'Incorrect (Fail)'};
boxplot([pass_corr_RT pass_incorr_RT fail_corr_RT fail_incorr_RT], group, 'Labels', labels);





















% 
% %% boxplot 생성
% f=figure;f.Position = [1645,857,829,594]; f.PaperSize = [21,30];
% hold on
% title('Response Time (for Correct and Incorrect Trials)')
% % xlabel('Subject Number')
% ylabel('RT(s)')
% 
% hc = boxplot([corr_RT; incorr_RT],[correct; incorrect], OutlierSize=10^(-200));
% idxc = find(ismember([SubInfoFile.Session], F));
% set(h(:,idxc),'Color','red');
% set(h(6,:),'Color','k');
% set(h(7,:),'MarkerEdgeColor','w');
% set(h(1:2,:),'LineStyle','-');
% 
% 
% boxplot([corr_RT(:) incorr_RT(:)], [repmat(1:size(corr_RT,2),1,size(corr_RT,1))';...
%     repmat(1:size(incorr_RT,2),1,size(incorr_RT,1))'+size(corr_RT,2)], 'labels',...
%     [cellstr(num2str((1:size(corr_RT,2))'))' cellstr(num2str((1:size(incorr_RT,2))'))']);
% 
% boxplot([corr_RT incorr_RT], 'labels', [cellstr(num2str((1:size(corr_RT,2))'))' cellstr(num2str((1:size(incorr_RT,2))'))']);
% 
% 
% % Remove NaN values from corr_RT
% corr_RT_noNaN = cell(size(corr_RT,2),1);
% for i = 1:size(corr_RT,2)
%     corr_RT_noNaN{i} = corr_RT(~isnan(corr_RT(:,i)),i);
% end
% % Create boxplot
% boxplot(cell2mat(corr_RT_noNaN'), 'labels', cellstr(num2str((1:size(corr_RT,2))')))
% title('RT for Correct Trials')
% ylabel('RT (s)')
% 
% 
% % Remove NaN values from incorr_RT
% incorr_RT_noNaN = cell(size(incorr_RT,2),1);
% for i = 1:size(incorr_RT,2)
%     incorr_RT_noNaN{i} = incorr_RT(~isnan(incorr_RT(:,i)),i);
% end
% 
% % Create boxplot
% boxplot([corr_RT_noNaN(:) incorr_RT_noNaN(:)], 'labels', {'Correct', 'Incorrect'})
% title('RT for Correct and Incorrect Trials')
% ylabel('RT (s)')
% 
% 
% 
% ax = gca;
% xTick = ax.XTick;
% xLim = ax.XLim;
% ylim([0 1.6])
% 
% for i = 1:length(idxc)
%     text(xTick(idxc(i)), yLim(1)-0.01*diff(yLim), ax.XTickLabel{idx(i)},...
%         'Color', 'red', 'HorizontalAlignment', 'center');
% ax.XTickLabel{idxc(i)} = '';
% end
% 


































