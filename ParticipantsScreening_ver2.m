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
allAccuracy = [];halfbias=[];allbias=[];halfAccuracy=[];
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    DataGroup.Correct.(varName) = Tc((Tc.Session(:) == Subs(i)),:);

%% Performance-Bias
Screening.(varName)= struct("Accuracy_all",[],"Accuracy_Lap",[],"Accuracy_Half",[],"RT",[],"Bias_all",[],"Bias_Lap",[],"Bias_Half",[]);
%Bias_overall trials
ButtonA = length(find(DataGroup.Overall.(varName).Choice==1));
ButtonB = length(find(DataGroup.Overall.(varName).Choice==2));
Screening.(varName).Bias_all = (ButtonA-ButtonB) / height(DataGroup.Overall.(varName));

%accuracy_overall trials
Screening.(varName).Accuracy_all = ((length(find(DataGroup.Correct.(varName).Correct==1)))/32)*100;

%Bias_per lap
Screening.(varName).Bias_Lap = struct;lapAccuracy =[];
for j=1:8
        lapnum = find(DataGroup.Overall.(varName).Lap==j);
        lapnumchoice = DataGroup.Overall.(varName).Choice(lapnum);

ButtonA_Lap = length(find(lapnumchoice==1));
ButtonB_Lap = length(find(lapnumchoice==2));
Screening.(varName).Bias_Lap(1,j) = (ButtonA_Lap-ButtonB_Lap) / height(lapnumchoice);

%accuracy_per lap
lapAccuracy = length(find(DataGroup.Correct.(varName).Correct(lapnum)==1));
Screening.(varName).Accuracy_Lap(1,j) = (lapAccuracy/4)*100;

%Bias_Last Half
if j==5
    halfchoice = DataGroup.Overall.(varName).Choice(lapnum(1):end);
    ButtonA_Half = length(find(halfchoice==1));
    ButtonB_Half = length(find(halfchoice==2));
    Screening.(varName).Bias_Half = (ButtonA_Half-ButtonB_Half) / height(halfchoice);

%accuracy_Last half
    halfcorr = DataGroup.Overall.(varName).Correct(lapnum(1):end);
    Screening.(varName).Accuracy_Half = (length(find(halfcorr))/ height(halfcorr))*100;

end
end

%% SubInfoFile에 정보 추가
allAccuracy = [allAccuracy; Screening.(varName).Accuracy_all];
halfAccuracy = [halfAccuracy; Screening.(varName).Accuracy_Half];
allbias = [allbias; Screening.(varName).Bias_all];
halfbias = [halfbias; Screening.(varName).Bias_Half];


%% 써먹지 않았지만 꽤 괜찮은 if문 (PASS/FAIL 구분)
% if ismember(Subs(i), P)
% elseif ismember(Subs(i), F)
% end


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
end;clear("temp")

%screening RT mean값 구해서 넣기
meanValues = nanmean(overall_RT)';
Screening.(varName).RT = meanValues(i);

end
SubInfoFile = addvars(SubInfoFile, allAccuracy, halfAccuracy, allbias, halfbias); 

%% RT boxplot 생성 - for each subject
figure('position',[1645 857 829 594]);
hold on
title('Response Time (for each Subject)')
xlabel('Subject')
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

%% RT boxplot 생성 - all sub, corr/incorr 2 boxes
figure('position',[1645 857 829 594]);
hold on
title('Response Time (for Correct and Incorrect Trials)')
xlabel('Correctness')
ylabel('RT(s)')

group = [ones(1,18), 2*ones(1,18)];
labels = {'Correct', 'Incorrect'};
xlim([0 3]);ylim([0 1.3]);
ha = boxplot([nanmean(corr_RT) nanmean(incorr_RT)],group, 'Labels', labels);


%% RT boxplot 생성 - PASS/Fail sub, corr/incorr 2boxes

figure('position',[1645 857 829 594]);
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

% boxplot
group = [ones(1,length(pass_corr_RT)), 2*ones(1,length(pass_incorr_RT)), ...
         3*ones(1,length(fail_corr_RT)), 4*ones(1,length(fail_incorr_RT))];
labels = {'Correct (Pass)', 'Incorrect (Pass)', 'Correct (Fail)', 'Incorrect (Fail)'};
boxplot([pass_corr_RT pass_incorr_RT fail_corr_RT fail_incorr_RT], group, 'Labels', labels);


%% RT boxplot 생성 - all sub, corr/incorr 2boxes

figure('position',[1300,700,1200,800]);

hold on
title('Response Time for Correct and Incorrect Trials by Subject')
xlabel('Subject')
ylabel('RT(s)')
set(h(1:2,:),'LineStyle','-');

% boxplot 그리기
group = [1:18, 19:36];
positions = [1:18, 19:36];
labels = {};
for i=1:18
    labels{end+1} = ['sub', num2str(SubInfoFile.Session(i)), '(Corr)'];
end
for i=1:18
    labels{end+1} = ['sub', num2str(SubInfoFile.Session(i)), '(Incorr)'];
end
boxplot([corr_RT incorr_RT], group, 'Labels', labels, 'Positions', positions, OutlierSize=10^(-200));

% 색상 변경하기
h = findobj(gca,'Tag','Box');
h = flipud(h);
for j=1:length(h)
    if group(j) <= 18
        patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.4);
        if find((j)==fail_group)
            lines = findobj(h(j), 'Type', 'Line');
        set(lines, 'Color', 'r', 'LineWidth', 2);
        end
    end
    if find(j-18==fail_group)
        lines = findobj(h(j), 'Type', 'Line');
        set(lines, 'Color', 'r', 'LineWidth', 2);
    end
end


























