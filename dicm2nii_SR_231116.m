clear all; clc;
addpath(genpath('\\147.47.203.154\LeeStorage1\E-Phys Analysis\fMRI_ocat\OCAT_DIR'));

%%
file_path_in = '../data/data_fmri_raw';
file_path_out = '../data/data_dicm2nii'; 

%% set subject's root directory
n_sbj = 31;


for sbj_i = 1: n_sbj
    c_sbj =  strcat("sub-",num2str(sbj_i, '%02.f'));
    
    %% path_in == raw data folder: T1 / field-mag / field-phase / func
    raw_sbj_folders = dir(fullfile(file_path_in, string(c_sbj), 'study*', 'series_*'));
    path_in = {};
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "T1"),1,"last"));
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"first"));
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"last"));
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "TASK"),1,"last"));
    %% make output directory
    path_sbj_out = fullfile(file_path_out,c_sbj);
     path_out = {};
     path_out{end+1} = fullfile(path_sbj_out,'anat');
     path_out{end+1} = fullfile(path_sbj_out,'fmap','mag');
     path_out{end+1} = fullfile(path_sbj_out,'fmap','phase');
     path_out{end+1} = fullfile(path_sbj_out,'func');
     
        if ~exist(path_sbj_out,"dir")
            mkdir(path_sbj_out);
            for i=1:4
                mkdir(path_out{i});
            end
        end
    
%% converting image to nifti
    for pth = 1:length(path_in)
    dicomDir = fullfile(path_in{1,pth}.folder, path_in{1,pth}.name);
    outputDir = path_out{1,pth};
   
   
    dicm2nii(dicomDir, outputDir, 'nii');
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
    'Name', 'Object-Context Association Task',...
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

