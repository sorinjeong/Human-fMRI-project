clear(sbj_events);
%% INPUT!!
n_sbj = 31; % enter the number of subjects

%% set event table
events=all_sbj_events;
sbj_info = sbj_info_file;
sbj_info = removevars(sbj_info,["Weight","Size"]);

% all_sbj_events=addvars(all_sbj_events,combi,repmat([1;2;3;4],(height(all_sbj_events)/4),1),repmat([1;2;1;2],(height(all_sbj_events)/4),1),'NewVariableNames',{'Combination','StopPoint','1324'});

sbj_perform=struct;
fn = {'all_accu','per_lap_accu','first_h_accu','second_h_accu','all_bias','per_lap_bias','first_h_bias','second_h_bias','rt_overall','rt_corr','rt_incorr','rt_mean'};
for i = 1:numel(fn)
    sbj_perform.(fn{i}) = [];
end
% all_accu=[]; per_lap_accu=[]; first_h_accu=[];second_h_accu=[];
% all_bias=[]; per_lap_bias=[]; first_h_bias=[];second_h_bias=[];

for sbj_i = 1: n_sbj
    c_sbj = strcat('sub', num2str(sbj_i, '%02.f'));
    disp(['Current subject: ', c_sbj]);

    % subject별 event 분리
    sbj_event.(c_sbj) = events((((sbj_i-1)*32)+1):sbj_i*32,:);

    %% Context-Object 조합
    combi_C = sbj_event.(c_sbj).Context_txt(sbj_event.(c_sbj).Association ==1);
    combi_O = sbj_event.(c_sbj).Obj_ID(sbj_event.(c_sbj).Association ==1);

    combi_FFCC = [sort(combi_O(find(combi_C=="F",2,"first"))); sort(combi_O(find(combi_C=="C",2,"first")))]';
    combi_FFCC = strjoin(string(combi_FFCC));
    disp(['combi_FFCC: ', combi_FFCC]);
    sbj_info.Combi_FFCC{sbj_i} = combi_FFCC;

%% Data Group (correct/overall)
data_group = struct("Correct",[],"Overall",[]);
data_group.Overall = sbj_event;
data_group.Correct = sbj_event;
data_group.Correct.(c_sbj)((data_group.Correct.(c_sbj).Correct_Num ~= 1),:) = [];

%% screening
for i = 1:numel(fn)
screening.(c_sbj).(fn{i}) = []; %struct('PASS',[],'all_accu',[],'per_lap_accu',[],'first_h_accu',[],'second_h_accu',[],'all_bias',[],'per_lap_bias',[],'first_h_bias',[],'second_h_bias',[]);
end
% accuracy_overall trials
screening.(c_sbj).all_accu = (height(data_group.Correct.(c_sbj))/32)*100;

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
screening.(c_sbj).rt_incorr = data_group.Overall.(c_sbj)((data_group.Overall.(c_sbj).Correct_Num == 2),'RT');
screening.(c_sbj).rt_mean = nanmean(screening.(c_sbj).rt_overall);


%% PASS/FAIL
if screening.(c_sbj).second_h_accu < 0.7 || screening.(c_sbj).second_h_bias > 0.2
    if screening.(c_sbj).second_h_accu < 0.7 && screening.(c_sbj).second_h_bias > 0.2
        sbj_info.PASS{sbj_i} = 1;
    else
        disp(['<second_half - ' c_sbj '>']), disp(['Accuracy: ', num2str(screening.(c_sbj).second_h_accu)]), disp(['Bias: ', num2str(screening.(c_sbj).second_h_bias)]);
        sbj_info.PASS{sbj_i} = input('Enter a value for PASS: ');
    end
else 
    sbj_info.PASS{sbj_i} = 0;
end

%% performance table 만들기
for i = 1:numel(fn)
        sbj_perform.(fn{i}) = [sbj_perform.(fn{i}); screening.(c_sbj).(fn{i})];
end






















