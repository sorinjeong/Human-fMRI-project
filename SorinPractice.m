% O-CAT data practice
%Input
NumOfSubs = 1;

%% subject numbering , folder root
for subname = 1:NumOfSubs
    Subjects{subname} = sprintf('Sub%.15g', (subname+84));
end
total_NumLogTable=[];
for fi = 1:numel(Subjects)
    Session=Subjects{fi};

% Responses -> 0 : incorrect trial + timeout, 1: correct trial
Responses = readmatrix(['Z:\E-Phys Analysis\fMRI_ocat\PilotData_analyzed\ver_230828\Total\LogTable_Total.xlsx'],'Sheet',Session,'Range','O2:O33');
Responses(Responses==2) = 0; Responses=Responses';

[pdata cdata] = LearningCurve_WinBugs_230829(Responses);

% Plot example
figure();
plotI(Responses, ones(1,length(Responses))); hold on;
title(['LearningCurve: ' Session],'FontSize',14,'FontWeight','bold')
xlabel('Trial Number');
ylabel('Pr(Correct Response)')
plot(pdata(:,2),'k-','LineWidth',1.2); hold on;
plot(pdata(:,3),'r-','LineWidth',1.2); hold on;
plot(pdata(:,4),'k-','LineWidth',1.2);
l=line([cdata cdata],[0 1]); l.Color='b'; l.LineWidth=1.8;




addpath(genpath('D:\leelab\Human fMRI projects\MATLAB\23.08.29_LearningCurve'))
saveas(gcf,['D:\leelab\Human fMRI projects\MATLAB\23.08.29_LearningCurve\Figures\' Session],'png');
hold off; close

end
