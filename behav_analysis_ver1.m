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
allAccuracy = [];halfbias_L=[];allbias=[];halfAccuracy_L=[];accu_perlap=[];bias_perlap=[];halfAccuracy_F=[];halfbias_F=[];
for i=1:length(Subs)
    varName = sprintf('SUB_%.15g',Subs(i));
    DataGroup.Correct.(varName) = Tc((Tc.Session(:) == Subs(i)),:);

%% Performance-Bias
Screening.(varName)= struct("Accuracy_all",[],"Accuracy_Lap",[],"Accuracy_Half_F",[],"Accuracy_Half_L",[],"RT",[],"Bias_all",[],"Bias_Lap",[],"Bias_Half_F",[],"Bias_Half_L",[]);
%Bias_overall trials
ButtonA = length(find(DataGroup.Overall.(varName).Choice==1));
ButtonB = length(find(DataGroup.Overall.(varName).Choice==2));
Screening.(varName).Bias_all = (ButtonA-ButtonB) / height(DataGroup.Overall.(varName));

%accuracy_overall trials
Screening.(varName).Accuracy_all = ((length(find(DataGroup.Correct.(varName).Correct==1)))/32)*100;

%Bias_per lap
lapAccuracy =[];
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
    halfchoice_L = DataGroup.Overall.(varName).Choice(lapnum(1):end);
    ButtonA_Half_L = length(find(halfchoice_L==1));
    ButtonB_Half_L = length(find(halfchoice_L==2));
    Screening.(varName).Bias_Half_L = (ButtonA_Half_L-ButtonB_Half_L) / height(halfchoice_L);

    halfchoice_F = DataGroup.Overall.(varName).Choice(1:lapnum(1)-1);
    ButtonA_Half_F = length(find(halfchoice_F==1));
    ButtonB_Half_F = length(find(halfchoice_F==2));
    Screening.(varName).Bias_Half_F = (ButtonA_Half_F-ButtonB_Half_F) / height(halfchoice_F);

%accuracy_Last half
    halfcorr_L = DataGroup.Overall.(varName).Correct(lapnum(1):end);
    Screening.(varName).Accuracy_Half_L = (length(find(rmmissing(halfcorr_L)))/ height(halfcorr_L))*100;

    halfcorr_F = DataGroup.Overall.(varName).Correct(1:lapnum(1)-1);
    Screening.(varName).Accuracy_Half_F = (length(find(rmmissing(halfcorr_F)))/ height(halfcorr_F))*100;
end
end

%% SubInfoFile에 정보 추가
allAccuracy = [allAccuracy; Screening.(varName).Accuracy_all];
halfAccuracy_L = [halfAccuracy_L; Screening.(varName).Accuracy_Half_L];
halfAccuracy_F = [halfAccuracy_F; Screening.(varName).Accuracy_Half_F];
allbias = [allbias; Screening.(varName).Bias_all];
halfbias_L = [halfbias_L; Screening.(varName).Bias_Half_L];
halfbias_F = [halfbias_F; Screening.(varName).Bias_Half_F];
accu_perlap = [accu_perlap Screening.(varName).Accuracy_Lap'];
bias_perlap = [bias_perlap Screening.(varName).Bias_Lap'];


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
SubInfoFile = addvars(SubInfoFile, allAccuracy,halfAccuracy_F, halfAccuracy_L, allbias, halfbias_F, halfbias_L); 

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%% accuracy plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fail_group = find(ismember([SubInfoFile.Session], F));
%% First + all + Last half accuracy % connecting the lines between subjects % first-half -> last-half
x=SubInfoFile.Session;
y_all=allAccuracy;
y_half_F=halfAccuracy_F;
y_half_L=halfAccuracy_L;

figure
hold on
title('Accuracy (for each Subject)',FontSize=14,FontWeight='bold')
xlabel('Subject')
ylabel('Accuracy (%)')

% Plot the data points for y_all, y_half and y_first
plot(x, y_all, 'k^', 'linewidth', 1, 'MarkerSize', 5);
plot(x, y_half_L,'k', 'marker',"diamond",'linestyle', 'none', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
plot(x, y_half_F,'k', 'marker',"diamond",'linestyle', 'none', 'MarkerSize', 5);

% Draw a line between each pair of data points
for i = 1:length(x)
    if find((i) == fail_group)
        line([x(i) x(i)], [y_half_F(i) y_half_L(i)], 'Color', 'r', 'LineStyle', '--')
    else
        line([x(i) x(i)], [y_half_F(i) y_half_L(i)], 'Color', 'k', 'LineStyle', '--')
    end
end

legend({'overall','last-half','first-half','','fail group'},'Location','southwest')


% fail group에 색칠하기
ax = gca;
ax.XTick = x;xTick=ax.XTick;ax.XLim = [84 103];
ax.XTickLabel = cellstr(num2str(x(:)));ylim([20 110]);yLim=ylim;

for i = 1:length(fail_group)
    if fail_group(i) <= length(xTick)
        text(xTick(fail_group(i)), yLim(1)-0.05*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center','Rotation',45);
        ax.XTickLabel{fail_group(i)} = '';
    end
end



%% accuracy for each lap
x = SubInfoFile.Session;
figure
hold on
title('Accuracy (for each Lap)',FontSize=14,FontWeight='bold')
xlabel('Subject')
ylabel('Accuracy (%)')
ylim([-10 110])

h=boxplot(accu_perlap,x, OutlierSize=1);
set(h(6,:),'Color','k','LineWidth',2);
set(h(1:2,:),'LineStyle','-');
% fail group에 색칠하기
ax = gca; xTick = ax.XTick; xLim = ax.XLim; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;
% 피험자번호
for i = 1:length(fail_group)
        text(xTick(fail_group(i)), yLim(1)-0.05*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center','Rotation',45);
    ax.XTickLabel{fail_group(i)} = '';
end
% box색깔
h = findobj(gca,'Tag','Box');
h = flipud(h);
for j=1:length(h)
    if find((j)==fail_group)
      lines = findobj(h(j), 'Type', 'Line');
      set(lines, 'Color', 'r', 'LineWidth', 1);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%% Bias plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First + all + Last bias % connect the dots for y_first-half and y_last-half
x = SubInfoFile.Session;
y_all = allbias;
y_half_L = halfbias_L;
y_half_F = halfbias_F;

figure
hold on
title('Bias Change (for each Subject)', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Subject')
ylabel('Button Press Bias')

% Plot the data points for y_all, y_half and y_first
plot(x, y_half_F, 'ko', 'linewidth', 1, 'MarkerSize', 5);
plot(x, y_half_L, 'ko', 'linewidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'k');

% Plot the data points for y_all with different marker based on the relationship between y_half_F and y_half_L
for i = 1:length(x)
    if y_half_F(i) > y_half_L(i)
        plot(x(i), y_all(i), 'kv', 'linewidth', 1, 'MarkerSize', 5,'HandleVisibility','off');
    else
        plot(x(i), y_all(i), 'k^', 'linewidth', 1, 'MarkerSize', 5,'HandleVisibility','off');
    end
end
plot(nan,nan,'k^');
% Draw a line between each pair of data points
for i = 1:length(x)
    if find((i) == fail_group)
        line([x(i) x(i)], [y_half_F(i) y_half_L(i)], 'Color', 'r', 'LineStyle', '--')
    else
        line([x(i) x(i)], [y_half_F(i) y_half_L(i)], 'Color', 'k', 'LineStyle', '--')
    end
end

legend({'first-half','last-half','overall','','fail group'},'Location','northeast')

% coloring fail group
ax = gca;
ax.XTick = x;
xTick = ax.XTick; ax.XLim = [84 103]; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;

for i = 1:length(fail_group)
    if fail_group(i) <= length(xTick)
        text(xTick(fail_group(i)), yLim(1)-0.05*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center','Rotation',45);
        ax.XTickLabel{fail_group(i)} = '';
    end
end

% Draw threshold at y=0
line(xlim,[0 0],'Color','k','LineStyle','--','HandleVisibility','off')

hold off


%% Horizontal!!! bias for each lap connecting the lines between subjects

% Create a new figure
x = SubInfoFile.Session;
figure('position',[1535 609 818 624]);
hold on
title('Bias (for each Lap)',FontSize=14,FontWeight='bold')
ylabel('Subject')
xlabel('Button Press Bias')

% Plot the horizontal box plot
h = boxplot(bias_perlap, 'Labels', x, 'orientation', 'horizontal');
set(h(6,:), 'Color', 'k', 'LineWidth', 2);
set(h(1:2,:),'LineStyle','-');

% Set the x-axis limits to be symmetric around zero
ax = gca;
xLim = max(abs(ax.XLim));
xlim([-xLim xLim]); ax.YTickLabel = cellstr(num2str(x(:)));
% coloring fail group

% 피험자번호
for i = 1:length(fail_group)
        text(ax.XLim(1)-0.02, ax.YTick(fail_group(i)), ax.YTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'right');
    ax.YTickLabel{fail_group(i)} = '';
end
% box색깔
h = findobj(gca,'Tag','Box');
h = flipud(h);
for j=1:length(h)
    if find((j)==fail_group)
      lines = findobj(h(j), 'Type', 'Line');
      set(lines, 'Color', 'r', 'LineWidth', 2);
    end
end




%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RT Plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RT boxplot 생성 - for each subject
figure('position',[1645 857 829 594]);
hold on
title('Response Time (for each Subject)',FontSize=14,FontWeight='bold')
xlabel('Subject')
ylabel('RT(s)')

h = boxplot(overall_RT,SubInfoFile.Session, OutlierSize=10^(-200));
fail_group = find(ismember([SubInfoFile.Session], F));
set(h(:,fail_group),'Color','red');
set(h(6,:),'Color','k');
set(h(7,:),'MarkerEdgeColor','w');
set(h(1:2,:),'LineStyle','-');

%fail group만 색칠하기!
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


%% Screening_RT plot - Corr/inCorr

%% RT boxplot 생성 - all sub, corr/incorr 2 boxes
figure('position',[1645 857 829 594]);
hold on
title('Response Time (for Correct and Incorrect Trials)',FontSize=14,FontWeight='bold')
xlabel('Correctness')
ylabel('RT(s)')

group = [ones(1,18), 2*ones(1,18)];
labels = {'Correct', 'Incorrect'};
xlim([0 3]);ylim([0 1.3]);
ha = boxplot([nanmean(corr_RT) nanmean(incorr_RT)],group, 'Labels', labels);


%% RT boxplot 생성 - PASS/Fail sub, corr/incorr 2boxes

figure('position',[1645 857 829 594]);
hold on
title('Response Time (for Correct and Incorrect Trials) by Group(P/F)',FontSize=14,FontWeight='bold')
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
title('Response Time for Correct and Incorrect Trials by Subject',FontSize=14,FontWeight='bold')
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


%%%%%%%%%%%%%%%%%%%%%%%%여기서부턴 FAIL group은 분석에서 제외! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Lap별로 피험자들의 RT 변화를 볼 수 있는 그래프
% Change in RT over Laps
RT_perLap=[];pass_RT=overall_RT; pass_RT(:,fail_group)=[];
for i=0:7; RT_perLap = [RT_perLap;nanmean((pass_RT((i*4)+1:(i+1)*4,:)))];end
RT_allSub_perLap = nanmean(RT_perLap,2);

figure
hold on
plot(RT_allSub_perLap,'k-o');
boxplot(RT_perLap')
title('Change in RT over Laps',FontSize=14,FontWeight='bold')
xlabel('Lap')
ylabel('RT(s)')

legend('Subject Average','Location','northeast')

hold off


%% Context별 accuracy와 RT를 lap별로, 피험자별로 확인할 수 있는 그래프 
% 1. Context - accuracy
Forest_accu=[];City_accu=[];Forest_RT=[];City_RT=[];
for i=1:length(pass_group)
    varName = sprintf('SUB_%.15g',Subs(pass_group(i)));
    FA = mean(DataGroup.PASS.(varName).Correct(find(DataGroup.PASS.(varName).Context==1)))*100;
    CA = mean(DataGroup.PASS.(varName).Correct(find(DataGroup.PASS.(varName).Context==2)))*100;
    FR = mean(DataGroup.PASS.(varName).RT(find(DataGroup.PASS.(varName).Context==1)));
    CR = mean(DataGroup.PASS.(varName).RT(find(DataGroup.PASS.(varName).Context==2)));    

Forest_accu = [Forest_accu FA];
City_accu = [City_accu CA];
Forest_RT = [Forest_RT FR];
City_RT = [City_RT CR];

end



% bidirectional plot
figure
hold on
subplot(2,1,1)
bh = barh([-City_accu' Forest_accu']);
set(bh(2),'FaceColor','g')
set(gca,'YDir','reverse')
title('Context-Dependent Accuracy Differences', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Subject')
xlabel('Accuracy(%)')
legend({'City','Forest'},'Location','bestoutside')


subplot(2,1,2)
bh = barh([-City_RT' Forest_RT']);
set(bh(2),'FaceColor','g')
set(gca,'YDir','reverse')
title('Context-Dependent RT Differences', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Subject')
xlabel('RT(s)')
legend({'City','Forest'},'Location','bestoutside')
hold off


% line plot
figure
subplot(2,1,1)
hold on
plot(Forest_accu,'g')
plot(City_accu,'b')
title('Context-Dependent Accuracy Differences', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Subject')
ylabel('Accuracy(%)')
legend({'Forest','City'},'Location','bestoutside')
xlim([0 13])

subplot(2,1,2)
hold on
plot(Forest_RT,'g')
plot(City_RT,'b')
hold off
title('Context-Dependent RT Differences', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Subject')
ylabel('RT(s)')
legend({'Accuracy','RT'},'Location','bestoutside')
xlim([0 13])


% Forest에서 city 뺀거, 그 차이를 넣어야하나??
figure
hold on
% yyaxis left
plot(Forest_accu - City_accu,'r')
ylabel('Accuracy(%)')

% yyaxis right
plot(Forest_RT - City_RT,'k')
ylabel('RT(s)')

title('Context-Dependent Accuracy Differences', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Subject')
legend({'Accuracy','RT'},'Location','bestoutside')
















