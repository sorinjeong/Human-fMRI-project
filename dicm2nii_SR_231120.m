clear all; clc;
addpath(genpath('\\147.47.203.154\LeeStorage1\E-Phys Analysis\fMRI_ocat\OCAT_DIR'));

%%
file_path_in = '../data/data_fmri_raw';
file_path_out = '../data/data_dicm2nii_1120'; 

%% set subject's root directory
n_sbj = 31;

%% new file name
new_prefixes = {'_T1w', '_magnitude1', '_magnitude2', '_phase', '_bold'};

for sbj_i = 1: n_sbj
    c_sbj =  strcat("sub-",num2str(sbj_i, '%02.f'));
    
    %% path_in == raw data folder: T1 / field-mag / field-phase / func
    raw_sbj_folders = dir(fullfile(file_path_in, string(c_sbj), 'study*', 'series_*'));
    path_in = {};
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "T1"),1,"last")); % T1 image
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"first")); % field map - magnitude
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"last")); % field map - phase
    path_in{end+1} = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "TASK"),1,"last")); % MR image
    %% make output directory
    path_sbj_out = fullfile(file_path_out,c_sbj,'ses-01');
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
    
    %파일 확장자별로 불러오기
    files_nii = dir(fullfile(outputDir, '*.nii'));
    files_json = dir(fullfile(outputDir, '*.json'));
    
        % 모든 변환된 파일들을 하나의 배열로 합칩니다.
    files = [files_nii; files_json];
    
        for i = 1:length(files)
            % 원본 파일 이름
            old_name = files(i).name;
        
            % 파일 확장자를 가져오기
            [~, ~, ext] = fileparts(old_name);
        
            % 새 파일 이름 생성
            if pth == 2
                if contains(old_name, 'e1')
                    new_name = [char(c_sbj), '_ses-01', '_task-OCAT', new_prefixes{2}, ext];
                else
                    new_name = [char(c_sbj), '_ses-01', '_task-OCAT', new_prefixes{3}, ext];
                end
            elseif pth == 1
                new_name = [char(c_sbj), '_ses-01', '_task-OCAT', new_prefixes{pth}, ext];
            else
                new_name = [char(c_sbj), '_ses-01', '_task-OCAT', new_prefixes{pth+1}, ext];
            end
            
            % 파일 이름 변경
            movefile(fullfile(outputDir, old_name), fullfile(outputDir, new_name));
        end
        
    % matfile 이름 변경
    files_mat = dir(fullfile(outputDir, '*.mat'));
    old_matfile = files_mat.name;
    new_matfile = [char(c_sbj), old_matfile];
    movefile(fullfile(outputDir, old_matfile), fullfile(outputDir, new_matfile));
        
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

