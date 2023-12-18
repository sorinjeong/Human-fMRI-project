% O-CAT data 231218
clear all; clc;

%% set path 
path_in= '../../data/data_learning_curve/Responses';
path_out='../../data/data_learning_curve';
addpath(genpath(path_out));

%Input
n_sbj = 31;
cdata_table = table('Size',[n_sbj 2], 'VariableTypes',["string" "double"],'VariableNames',["Session","acquisition_onset"]);



%% subject numbering , folder root
for sbj_i = 1: n_sbj
    c_sbj = strcat('sub-', num2str(sbj_i, '%02.f'));

% Responses -> 0 : incorrect trial + timeout, 1: correct trial
load(fullfile(path_in, [c_sbj '_Responses.mat']));

[pdata, cdata] = LearningCurve_WinBugs(Responses);
cdata_table.Session(sbj_i) = c_sbj;cdata_table.acquisition_onset(sbj_i) = cdata;

% Plot
figure();
plotI(Responses, ones(1,length(Responses))); hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')
plot(pdata(:,2),'k-','LineWidth',1.2); hold on;
plot(pdata(:,3),'r-','LineWidth',1.2); hold on;
plot(pdata(:,4),'k-','LineWidth',1.2); hold on;
ly=line([0 32], [0.5 0.5]);set(ly, 'LineStyle', '--', 'LineWidth', 1);hold on;
l=line([cdata cdata],[0 1]); l.Color='b'; l.LineWidth=1.8;

% Fill the area between plot(pdata(:,2)) and plot(pdata(:,3)) with a transparent gray color with 60% opacity
x = 1:length(pdata(:,2));
fill([x fliplr(x)], [pdata(:,2)' fliplr(pdata(:,4)')], [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

title(['LearningCurve: ' c_sbj],'FontSize',14,'FontWeight','bold');

% Remove the top and right axes lines
box off;


% save
saveas(gcf,[path_out '\' c_sbj '_learning_curve'],'png');
hold off; close

disp(['Completed processing for subject: ', c_sbj]);
end

writetable(cdata_table,[path_out '\AcquisitionOnset.xlsx']);

