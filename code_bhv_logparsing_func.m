clear all; clc; close all;
addpath(genpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV'));
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\code');


%% INPUT!!
task_name = 'OCAT'; %put the task name of your project
n_sbj = 31; % enter the number of subjects
is_save_output = 0; % if you want to save the output, type 1
is_open_plot = 0; % if you want to open the performance plot, type 1 / off -> No 4-group figure
create_dir = 0; % if you want to create the new directory (path), type 1 / default = 0
create_bhv_plots = 1; % if you want to perform behavior analysis , type 1
plot_save_output = 1; % if you want to save the behavior analysis result plot, type 1

%% set path
log_path_in = '../data/data_bhv_raw';
log_path_out = '../data/data_bhv_log_table'; 
plot_path_out = '../data/data_bhv_plot';
curve_path_out = '../data/data_learning_curve';
bids_path_in = 'D:\fMRI\OCAT_DIR\data\data_fmri_bids';
addpath(genpath(bids_path_in));


path = {log_path_in,log_path_out,plot_path_out,bids_path_in,curve_path_out};
if create_dir ==1
for i=1:length(path)
         mkdir(path{i});
end
end
%% subjects information file
sbj_info_file = readtable('../../OCAT_DIR/data/data_fmri_bids/participants.tsv','FileType','text');
sbj_info_file = removevars(sbj_info_file,["Weight","Size"]);



%% Start for loop
all_sbj_events = [];num_sbj_events=[];sbj_info_file_temp=sbj_info_file;
for sbj_i = 1: n_sbj
    c_sbj = strcat('sub-', num2str(sbj_i, '%02.f'));
    disp(['Current subject: ', c_sbj]);

 [all_sbj_events_temp,num_sbj_events_temp,fig,sbj_info_file_temp] = func_bhv_logparsing(path,task_name,sbj_i,is_save_output,is_open_plot,create_dir,sbj_info_file_temp);
   
disp(['Completed processing for subject: ', c_sbj]);
num_sbj_events = [num_sbj_events;num_sbj_events_temp];
all_sbj_events = [all_sbj_events;all_sbj_events_temp];
end
sbj_info_file=sbj_info_file_temp;

%%
combi= strcat(string(all_sbj_events.Context_txt),num2str(all_sbj_events.Obj_ID));
all_sbj_events=addvars(all_sbj_events,combi,repmat([1;2;3;4],(height(all_sbj_events)/4),1),repmat([1;2;1;2],(height(all_sbj_events)/4),1),'NewVariableNames',{'Combination','StopPoint','1324'});

if is_save_output == 1
writetable(all_sbj_events,fullfile(path{2},'total', 'all_sbj_events.csv'));
save(fullfile(path{2},'total','all_sbj_events') ,"all_sbj_events",'-mat');
writetable(num_sbj_events,fullfile(path{2},'total','num_sbj_events.csv'));
save(fullfile(path{2},'total','num_sbj_events') ,"num_sbj_events",'-mat');
writetable(sbj_info_file,fullfile(path{2},'total','sbj_info.xlsx'));
end

%% display messages
if is_save_output == 1
    disp('All tasks completed.');
    disp(['Outputs saved in: ', log_path_out]);
else
    disp('All tasks completed. No outputs were saved.');
end
disp(['Subjects processed: sub-01 to ', c_sbj]);


%% behavior analysis - performance plot!!
if create_bhv_plots ==1
    % Make output directory
    path_out = {};
    path_out{end+1} = fullfile(plot_path_out,'accuracy');
    path_out{end+1} = fullfile(plot_path_out,'bias');
    path_out{end+1} = fullfile(plot_path_out,'rt');

    if ~exist(path_out{3},"dir")
        for i=1:length(path_out)
            mkdir(path_out{i});
        end
    end

[sbj_perform, sbj_info, figs] = func_bhv_performance_plots(n_sbj, is_save_output, all_sbj_events, sbj_info_file);

    if plot_save_output == 1
       saveas(figs.plot1, fullfile(path_out{1}, 'Accuracy_for each subject_line.png'));
       saveas(figs.plot2, fullfile(path_out{1}, 'Accuracy_for each subject_box.png'));
       saveas(figs.plot3, fullfile(path_out{2}, 'Bias_for each subject_line.png'));
       saveas(figs.plot4, fullfile(path_out{2}, 'Bias_for each subject_box.png'));
       saveas(figs.plot5, fullfile(path_out{3}, 'RT_for each subject_line.png'));
       saveas(figs.plot6, fullfile(path_out{3}, 'RT_change_box_trial.png'));
       saveas(figs.plot7, fullfile(path_out{3}, 'RT_change_box_lap.png'));
       saveas(figs.plot8, fullfile(path_out{1}, 'Accuracy_change_line_trial.png'));
       saveas(figs.plot9, fullfile(path_out{1}, 'Accuracy_change_box_lap.png'));
       saveas(figs.plot10, fullfile(path_out{2}, 'Bias_change_box_lap.png'));    
       saveas(figs.plot11, fullfile(path_out{1}, 'Accuracy_change_line_Adjacent_trial.png'));        
    end
end