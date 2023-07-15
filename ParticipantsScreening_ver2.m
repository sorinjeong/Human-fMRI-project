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
allAccuracy = [];halfbias=[];allbias=[];halfAccuracy=[];accu_perlap=[];bias_perlap=[];
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
accu_perlap = [accu_perlap Screening.(varName).Accuracy_Lap'];
bias_perlap = [bias_perlap Screening.(varName).Bias_Lap'];

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

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%% accuracy plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% all + half accuracy
x=SubInfoFile.Session;
y_all=allAccuracy;
y_half=halfAccuracy;

figure
hold on
title('Accuracy (for each Subject)',FontSize=14,FontWeight='bold')
xlabel('Subject')
ylabel('Accuracy (%)')
plot(x,y_all,'k-o', x,y_half,'b--*','linewidth',1.5, 'MarkerSize',5)
% Plot the data points that are in fail_group
plot(x(fail_group),y_all(fail_group),'r*', x(fail_group),y_half(fail_group),'r*','MarkerSize',6);

legend({'overall','last-half','fail group'},'Location','southwest')


% fail group에 색칠하기
fail_group = find(ismember([SubInfoFile.Session], F));

ax = gca;
ax.XTick = x;
ax.XTickLabel = cellstr(num2str(x(:)));

for i = 1:length(fail_group)
    if fail_group(i) <= length(xTick)
        text(xTick(fail_group(i)), yLim(1)-0.05*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center','Rotation',45);
        ax.XTickLabel{fail_group(i)} = '';
    end
end


%% accuracy for each lap

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

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bias Plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Bidirectional Bar Chart
%% bias for each lap connecting the lines between subjects

x = SubInfoFile.Session;
y_all = allbias;
y_half = halfbias;

% Calculate the maximum value of the data
maxValue = max(max(abs([y_all(:); y_half(:)])));

% Create a new figure
figure('position',[1535 609 818 624]);
hold on
title('Bias (for each Subject)', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Button Press Bias')
ylabel('Subject')

% Plot the first data series as a horizontal bar chart
barh(x, y_all, 'FaceColor', [0.2 0.6 1])

% Plot the second data series as a horizontal bar chart, on top of the first data series
barh(x, y_half, 'FaceColor', [1 0.4 0.2], 'BarWidth', 0.5)

% Set the x-axis limits
xlim([-maxValue maxValue])
ylim([min(x)-1 max(x)+1])

% Add a legend
legend({'overall', 'last-half'}, 'Location', 'bestoutside')

% Add y-axis tick labels
set(gca, 'YTick', x)
set(gca, 'YTickLabel', cellstr(num2str(x(:))))

% Color fail group labels red
fail_group = find(ismember([SubInfoFile.Session], F));
for i = 1:length(fail_group)
    if fail_group(i) <= length(x)
        text(-maxValue-0.02*diff(xlim), x(fail_group(i)), num2str(x(fail_group(i))),...
            'Color', 'red', 'HorizontalAlignment', 'center');
         x(fail_group(i)) = '';
    end
end



%% bias for each lap connecting the lines between subjects

% Create a new figure
figure('position',[1535 609 818 624]);
hold on
title('Bias (for each Lap)',FontSize=14,FontWeight='bold')
ylabel('Subject')
xlabel('Button Press Bias')

% Plot the horizontal box plot
h=boxplot(bias_perlap',x,'orientation', 'horizontal');
set(h(6,:),'Color','k','LineWidth',2);
set(h(1:2,:),'LineStyle','-');

% Set the x-axis limits to be symmetric around zero
ax = gca;
xLim = max(abs(ax.XLim));
xlim([-xLim xLim]); ax.YTickLabel = cellstr(num2str(x(:)));
% coloring fail group

% 피험자번호
for i = 1:length(fail_group)
        text(ax.XLim(1), ax.YTick(fail_group(i)), ax.YTickLabel{fail_group(i)},...
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



%% bias for each lap connecting the lines between subjects

figure('position',[1535 609 818 624]);
hold on
title('Bias (for each Lap)',FontSize=14,FontWeight='bold')
xlabel('Subject')
ylabel('Button Press Bias')

h=boxplot(bias_perlap,x);
set(h(6,:),'Color','k','LineWidth',2);
set(h(1:2,:),'LineStyle','-');
% coloring fail group
ax = gca; xTick = ax.XTick; xLim = ax.XLim; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;
% 피험자번호
for i = 1:length(fail_group)
        text(xTick(fail_group(i)), yLim(1)-0.025*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center');
    ax.XTickLabel{fail_group(i)} = '';
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
































%%%%%% box plot

% %% all + half bias % connecting the lines between subjects
% x=SubInfoFile.Session;
% y_all=allbias;
% y_half=halfbias;
% 
% figure('position',[1535 609 818 624]);
% hold on
% title('Bias (for each Subject)',FontSize=14,FontWeight='bold')
% xlabel('Subject')
% ylabel('Button Press Bias')
% plot(x,y_all,'k-o', x,y_half,'b--*','linewidth',1.5, 'MarkerSize',5)
% % Plot the data points that are in fail_group
% plot(x(fail_group),y_all(fail_group),'r*', x(fail_group),y_half(fail_group),'r*','MarkerSize',6);
% 
% legend({'overall','last-half','fail group'},'Location','northeast')
% 
% 
% % coloring fail group
% % fail_group = find(ismember([SubInfoFile.Session], F));
% ax = gca;
% ax.XTick = x;
% xTick = ax.XTick; xLim = ax.XLim; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;
% 
% for i = 1:length(fail_group)
%     if fail_group(i) <= length(xTick)
%         text(xTick(fail_group(i)), yLim(1)-0.025*diff(yLim), ax.XTickLabel{fail_group(i)},...
%             'Color', 'red', 'HorizontalAlignment', 'center');
%         ax.XTickLabel{fail_group(i)} = '';
%     end
% end


%% all + half bias % connect the dots for y_all and y_half
x = SubInfoFile.Session;
y_all = allbias;
y_half = halfbias;

figure
hold on
title('Bias Change for each Subject', FontSize=14, FontWeight='bold')
xlabel('Subject')
ylabel('Button Press Bias')

% Plot the data points for y_all and y_half
plot(x,y_all,'ko', x,y_half,'bo','linewidth',1.5, 'MarkerSize',5, 'MarkerFaceColor','k');
plot(x,y_half,'bo','linewidth',1.5, 'MarkerSize',5, 'MarkerFaceColor','b');
% Draw a line between each pair of data points
for i = 1:length(x)
    if find((i)==fail_group)
        line([x(i) x(i)], [y_all(i) y_half(i)], 'Color', 'r', 'LineStyle', '--')
    else
        line([x(i) x(i)], [y_all(i) y_half(i)], 'Color', 'k', 'LineStyle', '--')
    end
end

% Create a dummy plot with a red dashed line
% p3=plot(NaN,NaN,'r--');

legend({'overall','last-half','','','fail group'},'Location','northeast')
% legend([p3],{'fail group'},'Location','northeast')

% coloring fail group
% fail_group = find(ismember([SubInfoFile.Session], F));
ax = gca;
ax.XTick = x;
xTick = ax.XTick; xLim = ax.XLim; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;

for i = 1:length(fail_group)
    if fail_group(i) <= length(xTick)
        text(xTick(fail_group(i)), yLim(1)-0.05*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center','Rotation',45);
        ax.XTickLabel{fail_group(i)} = '';
    end
end










%% bias for each lap connecting the lines between subjects

figure('position',[1535 609 818 624]);
hold on
title('Bias (for each Lap)',FontSize=14,FontWeight='bold')
xlabel('Subject')
ylabel('Button Press Bias')

h=boxplot(bias_perlap,x);
set(h(6,:),'Color','k','LineWidth',2);
set(h(1:2,:),'LineStyle','-');
% coloring fail group
ax = gca; xTick = ax.XTick; xLim = ax.XLim; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;
% 피험자번호
for i = 1:length(fail_group)
        text(xTick(fail_group(i)), yLim(1)-0.025*diff(yLim), ax.XTickLabel{fail_group(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center');
    ax.XTickLabel{fail_group(i)} = '';
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


% %% RT boxplot 생성 - for each subject
% figure('position',[1645 857 829 594]);
% hold on
% title('Response Time (for each Subject)',FontSize=14,FontWeight='bold')
% xlabel('Subject')
% ylabel('RT(s)')
% 
% h = boxplot(overall_RT,SubInfoFile.Session, OutlierSize=10^(-200));
% fail_group = find(ismember([SubInfoFile.Session], F));
% set(h(:,fail_group),'Color','red');
% set(h(6,:),'Color','k');
% set(h(7,:),'MarkerEdgeColor','w');
% set(h(1:2,:),'LineStyle','-');
% 
% %fail group만 색칠하기!
% ax = gca;
% xTick = ax.XTick;
% xLim = ax.XLim;
% ylim([0 1.6]);
% yLim = ylim;
% 
% for i = 1:length(fail_group)
%     text(xTick(fail_group(i)), yLim(1)-0.03*diff(yLim), ax.XTickLabel{fail_group(i)},...
%         'Color', 'red', 'HorizontalAlignment', 'center');
% ax.XTickLabel{fail_group(i)} = '';
% end
% 
% % % plot(meanValues) % mean값 넣고싶다면 이 주석을 풀어!
% % meanValues = nanmean(overall_RT);
% % for i = 1:18
% %     text(i, 1.55, jjnum2str(meanValues(i),2), 'HorizontalAlignment', 'center');
% % end
% 
% 
% %% Screening_RT plot - Corr/inCorr
% 
% %% RT boxplot 생성 - all sub, corr/incorr 2 boxes
% figure('position',[1645 857 829 594]);
% hold on
% title('Response Time (for Correct and Incorrect Trials)',FontSize=14,FontWeight='bold')
% xlabel('Correctness')]
% ylabel('RT(s)')
% 
% group = [ones(1,18), 2*ones(1,18)];
% labels = {'Correct', 'Incorrect'};
% xlim([0 3]);ylim([0 1.3]);
% ha = boxplot([nanmean(corr_RT) nanmean(incorr_RT)],group, 'Labels', labels);
% 
% 
% %% RT boxplot 생성 - PASS/Fail sub, corr/incorr 2boxes
% 
% figure('position',[1645 857 829 594]);
% hold on
% title('Response Time (for Correct and Incorrect Trials) by Group(P/F)',FontSize=14,FontWeight='bold')
% xlabel('Correctness')
% ylabel('RT(s)')
% xlim([0 5]);ylim([0 1.3]);
% pass_group = find(ismember([SubInfoFile.Session], P));
% 
% % 각 그룹의 데이터 계산하기
% pass_corr_RT = nanmean(corr_RT(:, pass_group));
% pass_incorr_RT = nanmean(incorr_RT(:, pass_group));
% fail_corr_RT = nanmean(corr_RT(:, fail_group));
% fail_incorr_RT = nanmean(incorr_RT(:, fail_group));
% 
% % boxplot
% group = [ones(1,length(pass_corr_RT)), 2*ones(1,length(pass_incorr_RT)), ...
%          3*ones(1,length(fail_corr_RT)), 4*ones(1,length(fail_incorr_RT))];
% labels = {'Correct (Pass)', 'Incorrect (Pass)', 'Correct (Fail)', 'Incorrect (Fail)'};
% boxplot([pass_corr_RT pass_incorr_RT fail_corr_RT fail_incorr_RT], group, 'Labels', labels);
% 
% 
% %% RT boxplot 생성 - all sub, corr/incorr 2boxes
% 
% figure('position',[1300,700,1200,800]);
% 
% hold on
% title('Response Time for Correct and Incorrect Trials by Subject',FontSize=14,FontWeight='bold')
% xlabel('Subject')
% ylabel('RT(s)')
% set(h(1:2,:),'LineStyle','-');
% 
% % boxplot 그리기
% group = [1:18, 19:36];
% positions = [1:18, 19:36];
% labels = {};
% for i=1:18
%     labels{end+1} = ['sub', num2str(SubInfoFile.Session(i)), '(Corr)'];
% end
% for i=1:18
%     labels{end+1} = ['sub', num2str(SubInfoFile.Session(i)), '(Incorr)'];
% end
% boxplot([corr_RT incorr_RT], group, 'Labels', labels, 'Positions', positions, OutlierSize=10^(-200));
% 
% % 색상 변경하기
% h = findobj(gca,'Tag','Box');
% h = flipud(h);
% for j=1:length(h)
%     if group(j) <= 18
%         patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.4);
%         if find((j)==fail_group)
%             lines = findobj(h(j), 'Type', 'Line');
%         set(lines, 'Color', 'r', 'LineWidth', 2);
%         end
%     end
%     if find(j-18==fail_group)
%         lines = findobj(h(j), 'Type', 'Line');
%         set(lines, 'Color', 'r', 'LineWidth', 2);
%     end
% end
% 

























