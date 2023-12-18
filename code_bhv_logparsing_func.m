clear all; clc; close all;
addpath(genpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV'));
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\code');

%% set path
log_path_in = '../data/data_bhv_raw';
log_path_out = '../data/data_bhv_log_table'; 
plot_path_out = '../data/data_bhv_plot';
curve_path_out = '../data/data_learning_curve';
derivatives_path_in = 'D:\fMRI\OCAT_DIR\data\data_fmri_bids\derivatives';
addpath(genpath(derivatives_path_in));

path = {log_path_in,log_path_out,plot_path_out,derivatives_path_in,curve_path_out};

%% INPUT!!
task_name = 'OCAT'; %put the task name of your project
n_sbj = 5; % enter the number of subjects
is_save_output = 1; % if you want to save the output, type 1
is_open_plot = 1; % if you want to open the performance plot, type 1 / off -> No 4-group figure

%% Start for loop
all_sbj_events = [];num_sbj_events=[];
for sbj_i = 1: n_sbj
    c_sbj = strcat('sub-', num2str(sbj_i, '%02.f'));
    disp(['Current subject: ', c_sbj]);

 [all_sbj_events_temp,num_sbj_events_temp,fig] = func_bhv_logparsing(path,task_name,sbj_i,is_save_output,is_open_plot);
   
disp(['Completed processing for subject: ', c_sbj]);
num_sbj_events = [num_sbj_events;num_sbj_events_temp];
all_sbj_events = [all_sbj_events;all_sbj_events_temp];

end

%%
% Calculate the number of subplots needed
numFigs = numel(fig);
numSubplots = ceil(numFigs / 4);

% Loop over each subplot
for i = 1:numSubplots
    % Create a new figure for the subplot
    f = figure('Name', sprintf('sub-%02dto%02d_Performance.png', (i-1)*4+1, min(i*4, numFigs)));
    
    % Loop over each figure in the current subplot
    for j = 1:4
        % Calculate the current figure index
        figIdx = (i-1)*4 + j;
        
        % If the current figure index is less than or equal to the total number of figures
        if figIdx <= numFigs
            % Copy the contents of the current figure to the subplot
            s = subplot(2, 2, j, 'Parent', f);
            ax = findall(fig{figIdx}, 'type', 'axes');
            copyobj(allchild(ax), s);
        end
    end
    
    % Save the subplot figure
    saveas(f, fullfile(path{3}, get(f, 'Name')));
end


%%
combi= strcat(string(all_sbj_events.Context_txt),num2str(all_sbj_events.Obj_ID));
all_sbj_events=addvars(all_sbj_events,combi,repmat([1;2;3;4],(height(all_sbj_events)/4),1),repmat([1;2;1;2],(height(all_sbj_events)/4),1),'NewVariableNames',{'Combination','StopPoint','1324'});

if is_save_output == 1
writetable(all_sbj_events,fullfile(path{2},'total', 'all_sbj_events.csv'));
save(fullfile(path{2},'total','all_sbj_events') ,"all_sbj_events",'-mat');
writetable(num_sbj_events,fullfile(path{2},'total','num_sbj_events.csv'));
save(fullfile(path{2},'total','num_sbj_events') ,"num_sbj_events",'-mat');
end

%% display messages
if is_save_output == 1
    disp('All tasks completed.');
    disp(['Outputs saved in: ', log_path_out]);
else
    disp('All tasks completed. No outputs were saved.');
end
disp(['Subjects processed: sub-01 to ', c_sbj]);
