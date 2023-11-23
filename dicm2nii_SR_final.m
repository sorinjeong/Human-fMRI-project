clear all; clc;
addpath(genpath('\\147.47.203.154\LeeStorage1\E-Phys Analysis\fMRI_ocat\OCAT_DIR'));

%%
file_path_in = '../data/data_fmri_raw_Num';
file_path_out = '../data/data_fmri_bids'; 

%% set subject's root directory
n_sbj = 18;

for sbj_i = 2: n_sbj
    c_sbj = num2str(sbj_i, '%02.f');
    
    %% path_in == raw data folder: T1 / field-mag / field-phase / func
    raw_sbj_folders = dir(fullfile(file_path_in, c_sbj, 'study*', 'series_*'));
    path_in = {};
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "T1"),1,"last")); % T1 image
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "TASK"),1,"last")); % MR image
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"first")); % field map - magnitude
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"last")); % field map - phase

    % %% make output directory
    % path_sbj_out = fullfile(file_path_out,c_sbj,'ses-01');
    %  path_out = {};
    %  path_out{end+1} = fullfile(path_sbj_out,'anat');
    %  path_out{end+1} = fullfile(path_sbj_out,'func');
    %  path_out{end+1} = fullfile(path_sbj_out,'fmap');
    %  path_out{end+1} = fullfile(path_sbj_out,'fmap');
    % 
    %     if ~exist(path_sbj_out,"dir")
    %         mkdir(path_sbj_out);
    %         for i=1:3
    %             mkdir(path_out{i});
    %         end
    %     end
    % 
%% converting image to nifti
    for pth = 1:length(path_in)
    dicomDir = fullfile(path_in{1,pth}.folder, path_in{1,pth}.name);
    outputDir = file_path_out;  
    
   
    dicm2nii_OCAT(dicomDir, outputDir, 'bids');
    
    end
end


%     %% arrayfun / cellfun 사용!
%     for pth = 1:length(path_in)
%     file_list = dir(fullfile(path_in{1,pth}.folder, path_in{1,pth}.name));
%     file_name = arrayfun(@(x) x.name, file_list, 'uni',0);
%     flag_dicom = cellfun(@(x) contains(x, '.IMA'), file_name);
% 
%     file_list = arrayfun(@(x) fullfile(x.folder, x.name), file_list(flag_dicom), 'uni', 0);
%     

%% JSON file
% Make 'dataset_description.json' as structure
dataset_description = struct(...
    'Name', 'Object-Context Association Task','TaskName', 'OCAT',...
    'BIDSVersion', '1.0.2',...
    'License', 'CC0',...
    'Authors', {{'Sorin Jeong', 'Zyeon Kim'}},...
    'Acknowledgements', 'We extend our thanks to Ji-Sun Kim and Jae-Min Seol for their assistance in data acquisition and analysis. We also want to express our gratitude to all the participants of this study.',...
    'HowToAcknowledge', 'Please cite this paper.');

% save as JSON file
jsonStr = jsonencode(dataset_description);
fid = fopen([file_path_out '\dataset_description.json'], 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);

% ,...
%     'Funding', {{'Grant 1', 'Grant 2'}},...
%     'ReferencesAndLinks', {{'Link 1', 'Link 2'}},...
%     'DatasetDOI', '10.0.2.3/abcd'...


% 'Participant_id' 열에서 첫 번째 행을 제외하고 "sub-" 추가
cd(outputDir)
T = readtable('participants.tsv', 'FileType', 'text', 'Delimiter', '\t');
T.participant_id(2:end) = strcat('sub-', string(T.participant_id(2:end)));
writetable(T, 'participants.tsv', 'FileType', 'text', 'Delimiter', '\t');
