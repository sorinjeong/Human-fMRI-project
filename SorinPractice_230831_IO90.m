% O-CAT data practice

%Input
NumOfSubs = 18;
cdata_table = zeros(NumOfSubs,1);

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

[pdata cdata] = LearningCurve_WinBugs_230829_IO90(Responses);
cdata_table(fi) = cdata;

% Plot
figure();
plotI(Responses, ones(1,length(Responses))); hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')
plot(pdata(:,2),'k-','LineWidth',1.2); hold on;
plot(pdata(:,3),'r-','LineWidth',1.2); hold on;
plot(pdata(:,4),'k-','LineWidth',1.2);
l=line([cdata cdata],[0 1]); l.Color='b'; l.LineWidth=1.8;

% Fill the area between plot(pdata(:,2)) and plot(pdata(:,3)) with a transparent gray color with 60% opacity
x = 1:length(pdata(:,2));
fill([x fliplr(x)], [pdata(:,2)' fliplr(pdata(:,4)')], [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

% Specify the path to the .mat file
addpath(genpath('D:\leelab\Human fMRI projects\MATLAB\23.08.29_LearningCurve'))
mat_file_path = 'D:\leelab\Human fMRI projects\MATLAB\23.08.29_LearningCurve\P';

% Load the contents of the .mat file into a structure
data = load(mat_file_path);

% Access the contents of the structure
P = data.P;

% Add the subtitle
if ismember((fi+84), P)
    subtitle = 'PASS';
else
    subtitle = 'FAIL';
end
t = title(['LearningCurve: ' Session ' (' subtitle ')'],'FontSize',14,'FontWeight','bold');

% Remove the top and right axes lines
box off;


%save

saveas(gcf,['D:\leelab\Human fMRI projects\MATLAB\23.08.29_LearningCurve\Figures_IO90\' Session],'png');
hold off; close

end

