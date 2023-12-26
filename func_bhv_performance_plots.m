function [sbj_perform, sbj_info, figs] = func_bhv_performance_plots(n_sbj, is_save_output, all_sbj_events, sbj_info_file)

    % Set event table
    events=all_sbj_events;
    sbj_info = sbj_info_file;
    sbj_info.participant_id = cellfun(@(x) strrep(x, '-', ''), sbj_info.participant_id, 'UniformOutput', false);

    screening=struct; sbj_perform=struct;
    fn = {'all_accu','per_lap_accu','first_h_accu','second_h_accu','all_bias','per_lap_bias','first_h_bias','second_h_bias','rt_overall','rt_corr','rt_incorr','rt_mean'};
    fn_single = {'all_accu','first_h_accu','second_h_accu','all_bias','first_h_bias','second_h_bias','rt_mean'};
    fn_par = {'per_lap_accu','per_lap_bias'};
    overall_RT=[];

    for i = 1:numel(fn_single); sbj_perform.(fn_single{i}) = []; end
    for i = 1:numel(fn_par); box.(fn_par{i}) = []; end

    data_group = struct("Correct",[],"Overall",[]);
    for sbj_i = 1: n_sbj
        c_sbj = strcat('sub', num2str(sbj_i, '%02.f'));
        disp(['Current subject: ', c_sbj]);

        % Data Group (correct/overall)
        data_group.Overall = sbj_event;
        data_group.Correct = sbj_event;
        data_group.Correct.(c_sbj)((data_group.Correct.(c_sbj).Correct_Num ~= 1),:) = [];
        % screening
        for i = 1:numel(fn)
            screening.(c_sbj).(fn{i}) = []; %struct('PASS',[],'all_accu',[],'per_lap_accu',[],'first_h_accu',[],'second_h_accu',[],'all_bias',[],'per_lap_bias',[],'first_h_bias',[],'second_h_bias',[]);
        end
        % accuracy_overall trials
        screening.(c_sbj).all_accu = (height(data_group.Correct.(c_sbj))/32);

        % bias_overall_trials 
        choice = [data_group.Overall.(c_sbj).Choice_Num{:}]';
        button_A = length(find(choice==1));
        button_B = length(find(choice==2));
        screening.(c_sbj).all_bias(end+1) = (button_A-button_B)/32;

        for lap = 1:8
            % accuracy_per_lap
            screening.(c_sbj).per_lap_accu(end+1) = (length(find(data_group.Correct.(c_sbj).Lap == lap))/4);

            % bias_per_lap
            lap_idx = find(data_group.Overall.(c_sbj).Lap == lap);
            choice = [data_group.Overall.(c_sbj).Choice_Num{lap_idx}]';
            button_A = length(find(choice==1));
            button_B = length(find(choice==2));
            screening.(c_sbj).per_lap_bias(end+1) = (button_A-button_B)/4;
        end

        % accuracy_half
        screening.(c_sbj).first_h_accu = sum(screening.(c_sbj).per_lap_accu(1:4))/4;
        screening.(c_sbj).second_h_accu = sum(screening.(c_sbj).per_lap_accu(5:end))/4;

        % bias_half
        screening.(c_sbj).first_h_bias = sum(screening.(c_sbj).per_lap_bias(1:4))/4;
        screening.(c_sbj).second_h_bias = sum(screening.(c_sbj).per_lap_bias(5:end))/4;

        %% RT
        screening.(c_sbj).rt_overall = data_group.Overall.(c_sbj).RT(data_group.Overall.(c_sbj).Correct_Num ~= 1);
        screening.(c_sbj).rt_corr = data_group.Correct.(c_sbj).RT;
        screening.(c_sbj).rt_incorr = data_group.Overall.(c_sbj).RT(data_group.Overall.(c_sbj).Correct_Num == 2);
        screening.(c_sbj).rt_mean = nanmean(screening.(c_sbj).rt_overall);

        %% PASS/FAIL
        % 0 = fail / 1 = pass / input = need to decide
        if screening.(c_sbj).second_h_accu < 0.7 || screening.(c_sbj).second_h_bias > 0.2
            if screening.(c_sbj).second_h_accu < 0.7 && screening.(c_sbj).second_h_bias > 0.2
                sbj_info.PASS{sbj_i} = 0;
            else
                disp(['<second_half - ' c_sbj '>']), disp(['Accuracy: ', num2str(screening.(c_sbj).second_h_accu)]), disp(['Bias: ', num2str(screening.(c_sbj).second_h_bias)]);
                sbj_info.PASS{sbj_i} = input('Enter a value for PASS: ');
            end
        else 
            sbj_info.PASS{sbj_i} = 1;
        end

        %% performance table 만들기
        for i = 1:numel(fn_single)
            sbj_perform.(fn_single{i})(sbj_i,1) = screening.(c_sbj).(fn_single{i});
        end
        for i = 1:numel(fn_par)
            box.(fn_par{i})(:,sbj_i) = [screening.(c_sbj).(fn_par{i})];
        end
        overall_RT(:,sbj_i)=[data_group.Overall.(c_sbj).RT];
    end
    sbj_info.PASS=cell2mat(sbj_info.PASS);
    sbj_perform = struct2table(sbj_perform);
    sbj_perform = horzcat(sbj_info,sbj_perform);
    
    % save sbj_info
    if is_save_output == 1
        writetable(sbj_info,[plot_path_out '\sbj_info.csv']);
    end

    % Create figures
    figs = struct();

    % Plotting code goes here...


%% %%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%
fail_idx = find(sbj_perform.PASS == 0);
pass_idx = find(sbj_perform.PASS == 1);

close all;
%% First + all + Last half accuracy % connecting the lines between subjects % first-half -> last-half

x=[1:31]';
y_all=sbj_perform.all_accu;
y_half_f=sbj_perform.first_h_accu;
y_half_s=sbj_perform.second_h_accu;

figs.plot1 = figure(Position=[1000,520,1391,878]);
hold on
title('Accuracy (for each subject)','FontSize', 14, 'FontWeight', 'bold')
subtitle('overall / first- / second- half로 나누어 계산')
xlabel('Subject')
ylabel('Accuracy')




% Plot the data points for y_all, y_half and y_first
plot(x, y_half_f,'k', 'marker',"diamond",'linestyle', 'none', 'MarkerSize', 5);
plot(x, y_half_s,'k', 'marker',"diamond",'linestyle', 'none', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

% Plot the data points for y_all with different marker based on the relationship between y_half_f and y_half_s
for i = 1:length(x)
    if y_half_f(i) > y_half_s(i)
        plot(x(i), y_all(i), 'kv', 'linewidth', 1, 'MarkerSize', 6, 'HandleVisibility', 'off');
    else
        plot(x(i), y_all(i), 'k^', 'linewidth', 1, 'MarkerSize', 5, 'HandleVisibility', 'off');
    end
end
plot(nan, nan, 'k^');

% Draw a line between each pair of data points
for i = 1:length(x)
    if ismember(i, fail_idx)
        line([x(i) x(i)], [y_half_f(i) y_half_s(i)], 'Color', 'r', 'LineStyle', '--')
    else
        line([x(i) x(i)], [y_half_f(i) y_half_s(i)], 'Color', 'k', 'LineStyle', '--')
    end
end


legend({'first-half','last-half','overall','','fail group'},'Location','southwest')


% fail group에 색칠하기
ax = gca;
ax.XTick = x;
ax.XTickLabel = cellstr(num2str((1:31)'));
xlim([0 32]);
ylim([0 1]);yLim = ylim;


for i = 1:length(fail_idx)
        text(ax.XTick(fail_idx(i)), ax.YLim(1)-0.02, ax.XTickLabel{fail_idx(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center');
    ax.XTickLabel{fail_idx(i)} = '';
end

box on; hold off


%% accuracy for each lap

figs.plot2=figure(Position=[1000,520,1391,878]);
hold on
title('Accuracy (for each subject)',FontSize=14,FontWeight='bold')
subtitle('Lap별로 나누어 계산')
xlabel('Subject')
ylabel('Accuracy')

h=boxplot(box.per_lap_accu,x, OutlierSize=1);
set(h(6,:),'Color','k','LineWidth',2);
set(h(1:2,:),'LineStyle','-');
% fail group에 색칠하기
ax = gca; xTick = ax.XTick; xLim = ax.XLim; ax.XTickLabel = cellstr(num2str(x(:)));ylim([0 1]);yLim = ylim;
for i = 1:length(fail_idx)
        text(ax.XTick(fail_idx(i)), ax.YLim(1)-0.02, ax.XTickLabel{fail_idx(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center');
    ax.XTickLabel{fail_idx(i)} = '';
end

% box색깔
h = findobj(gca,'Tag','Box');
h = flipud(h);

      lines = findobj(h(fail_idx), 'Type', 'Line');
      set(lines, 'Color', 'r', 'LineWidth', 1);

box on; hold off




%% %%%%%%%%%%%%%%%%%%%%%%%% Bias plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First + all + Last bias % connect the dots for y_first-half and y_last-half
y_all=sbj_perform.all_bias;
y_half_f=sbj_perform.first_h_bias;
y_half_s=sbj_perform.second_h_bias;

figs.plot3=figure('position',[399,393,1072,839]);
hold on
title('Bias (for each subject)',FontSize=14,FontWeight='bold')
subtitle('overall / first- / second- half로 나누어 계산')
xlabel('Subject')
ylabel('Button-Pressing Bias')

% Plot the data points for y_all, y_half and y_first
plot(x, y_half_f, 'ko', 'linewidth', 1, 'MarkerSize', 5);
plot(x, y_half_s, 'ko', 'linewidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'k');

% Plot the data points for y_all with different marker based on the relationship between y_half_f and y_half_s
for i = 1:length(x)
    if y_half_f(i) > y_half_s(i)
        plot(x(i), y_all(i), 'kv', 'linewidth', 1, 'MarkerSize', 5,'HandleVisibility','off');
    else
        plot(x(i), y_all(i), 'k^', 'linewidth', 1, 'MarkerSize', 5,'HandleVisibility','off');
    end
end
plot(nan,nan,'k^');
% Draw a line between each pair of data points
for i = 1:length(x)
    if ismember(i, fail_idx)
        line([x(i) x(i)], [y_half_f(i) y_half_s(i)], 'Color', 'r', 'LineStyle', '--')
    else
        line([x(i) x(i)], [y_half_f(i) y_half_s(i)], 'Color', 'k', 'LineStyle', '--')
    end
end

legend({'first-half','last-half','overall','','fail group'},'Location','northeast')

% coloring fail group
ax = gca;
ax.XTick = x;
xTick = ax.XTick; ax.XLim = [0 32]; ax.XTickLabel = cellstr(num2str(x(:)));yLim = ylim;

for i = 1:length(fail_idx)
        text(ax.XTick(fail_idx(i)), ax.YLim(1)-0.02, ax.XTickLabel{fail_idx(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center');
    ax.XTickLabel{fail_idx(i)} = '';
end

% Draw threshold at y=0
line(xlim,[0 0],'Color','k','LineStyle','--','HandleVisibility','off')

box on; hold off

%% Horizontal!!! bias for each lap connecting the lines between subjects

% Create a new figure
figs.plot4 = figure('position',[571,137,820,1025]);
hold on
title('Bias (for each subject)',FontSize=14,FontWeight='bold')
subtitle('Lap별로 나누어 계산')
ylabel('Subject')
xlabel('Button-Pressing Bias')

% Plot the horizontal box plot
h = boxplot(box.per_lap_bias, 'Labels', x, 'orientation', 'horizontal');
set(h(6,:), 'Color', 'k', 'LineWidth', 2);
set(h(1:2,:),'LineStyle','-');

% Set the x-axis limits to be symmetric around zero
ax = gca;
xLim = max(abs(ax.XLim));
xlim([-xLim xLim]); ax.YTickLabel = cellstr(num2str(x(:)));
% coloring fail group

% 피험자번호
for i = 1:length(fail_idx)
        text(ax.XLim(1)-0.02, ax.YTick(fail_idx(i)), ax.YTickLabel{fail_idx(i)},...
            'Color', 'red', 'HorizontalAlignment', 'right');
    ax.YTickLabel{fail_idx(i)} = '';
end
% box색깔
h = findobj(gca,'Tag','Box');
h = flipud(h);
lines = findobj(h(fail_idx), 'Type', 'Line');
set(lines, 'Color', 'r', 'LineWidth', 2);


box on; hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% RT Plot!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT boxplot 생성 - for each Lap
figs.plot5=figure('position',[1645 857 829 594]);
hold on
title('Response Time (for each subject)',FontSize=14,FontWeight='bold')
subtitle('Lap별로 나누어 계산')
xlabel('Subject')
ylabel('RT(s)')

h = boxplot(overall_RT,x, OutlierSize=10^(-200));
set(h(:,fail_idx),'Color','red');
set(h(6,:),'Color','k');
set(h(7,:),'MarkerEdgeColor','w');
set(h(1:2,:),'LineStyle','-');

%fail group만 색칠하기!
ax = gca;
xTick = ax.XTick;
xLim = ax.XLim;
ylim([0 1.6]);
yLim = ylim;

for i = 1:length(fail_idx)
        text(ax.XTick(fail_idx(i)), ax.YLim(1)-0.02, ax.XTickLabel{fail_idx(i)},...
            'Color', 'red', 'HorizontalAlignment', 'center');
    ax.XTickLabel{fail_idx(i)} = '';
end


%% %%%%%%%%%%%%%%%%%%%%%%여기서부턴 FAIL group은 분석에서 제외! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 피험자들의 RT 변화를 볼 수 있는 그래프
% Change in RT over Trials
pass_RT=overall_RT(:,pass_idx);
figs.plot6=figure('Position',[874 447 1685 951]);
hold on
plot(nanmean(pass_RT,2),'k-o', LineWidth=2);
boxplot(pass_RT')
title('Change in RT over Trials',FontSize=14,FontWeight='bold')
xlabel('Trial')
ylabel('RT(s)')
xlim([0 33])
legend('Subject Average','Location','northeast')

% Add vertical lines
for i = 1:7
    ia=(i*4)+0.5;
    line([ia ia], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '-','HandleVisibility','off');
end

box on; hold off

%% 
% Change in RT over Laps **230828 예쁘게!
RT_perLap=[];
for i=0:7; RT_perLap = [RT_perLap;nanmean((pass_RT((i*4)+1:(i+1)*4,:)))];end
RT_pass_perLap = nanmean(RT_perLap,2);

figs.plot7=figure;
hold on
plot(RT_pass_perLap,'Color','#013C58','Marker','o','MarkerFaceColor','#013C58','LineWidth',1.6);
h = boxplot(RT_perLap');
set(h,{'linew'},{2})
set(h,{'Color'},{[0.0039 0.2353 0.3451]})

% Change the transparency of the box plot
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.0039 0.2353 0.3451],'FaceAlpha',0.5);
end

% Change the line width of the outliers
h = findobj(gca,'Tag','Outliers');
set(h,{'MarkerSize'},{2})

set(gca, 'box', 'off')


title('Change in RT over Laps','FontSize',14,'FontWeight','bold')
xlabel('Lap')
ylabel('RT(s)')

legend('Subject Average','Location','northeast')

box on; hold off



%% 피험자들의 Accuracy 변화를 볼 수 있는 그래프
% Change in Accuracy over Trials

overall_accu=[];
for sbj_i = 1: n_sbj
    c_sbj = strcat('sub', num2str(sbj_i, '%02.f'));
    temp = sbj_event.(c_sbj).Correct_Num;
    temp(temp== 2) = 0;
    overall_accu= [overall_accu temp];
end
pass_Accuracy=overall_accu(:,pass_idx);
figs.plot8=figure('Position',[874 447 1685 951]);
hold on
plot(nanmean(pass_Accuracy,2),'k-o', LineWidth=2);
% boxplot(pass_Accuracy')
title('Change in Accuracy over Trials',FontSize=14,FontWeight='bold')
xlabel('Trial')
ylabel('Accuracy')
xlim([0 33]); ylim([0 1])
legend('Subject Average','Location','southeast')

% Add vertical lines
for i = 1:7
    ia=(i*4)+0.5;
    line([ia ia], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '-','HandleVisibility','off');
end

box on; hold off


%%
% Change in accuracy over Laps _ **230828 더 예쁘게 
pass_accu_lap = box.per_lap_accu(:,pass_idx);
accu_pass_perLap=nanmean(pass_accu_lap,2);

figs.plot9=figure;
hold on
plot(accu_pass_perLap,'Color','#F5564E','Marker','o','MarkerFaceColor','#F5564E','LineWidth',1.6);
h = boxplot(pass_accu_lap');
set(h,{'linew'},{1.6})
set(h,{'Color'},{[0.9608 0.3373 0.3059]})
ylim([0 1])

% Change the transparency of the box plot
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.9608 0.3373 0.3059],'FaceAlpha',0.5);
end

% Change the color and line width of the outliers
h = findobj(gca,'Tag','Outliers');
set(h,{'MarkerSize'},{1.6})
set(h,{'MarkerEdgeColor'},{[0.2314 0.2471 0.2745]})

% Change the color and line width of the median line
h = findobj(gca,'Tag','Median');
set(h,{'Color'},{[0.396, 0.263, 0.129]})
set(h,{'LineWidth'},{3})

set(gca, 'box', 'off')


title('Change in Accuracy over Laps','FontSize',14,'FontWeight','bold')
xlabel('Lap')
ylabel('Accuracy')

legend('Subject Average','Location','southeast')

box on; hold off


%% Lap별로 피험자들의 Bias 변화를 볼 수 있는 그래프
% Change in bias over Laps

pass_bias_lap = abs(box.per_lap_bias(:,pass_idx));
bias_pass_perLap=nanmean(pass_bias_lap,2);

figs.plot10=figure;
hold on
plot(bias_pass_perLap,'k-o');
boxplot(pass_bias_lap')
title('Change in Bias over Laps',FontSize=14,FontWeight='bold')
xlabel('Lap')
ylabel('Bias')
ylim([-0.050 1.05])
legend('Subject Average','Location','northeast')

box on; hold off

    return
end
