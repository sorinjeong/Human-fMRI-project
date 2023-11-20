clear all; clc;
restoredefaultpath; % Restore the MATLAB path to its default state.
addpath('C:\Users\Leelab\Documents\MATLAB\spm12');
addpath('C:\Users\Leelab\Documents\MATLAB\spm12\external\fieldtrip'); % Add the FieldTrip directory back to the MATLAB path.
ft_defaults; % Reset the FieldTrip defaults.

addpath(genpath('\\147.47.203.154\LeeStorage1\E-Phys Analysis\fMRI_ocat\OCAT_DIR'));

%%
spm('defaults','FMRI');
spm_jobman('initcfg');

%%
file_path_in = '../data/data_fmri_raw';
file_path_out = '../data/data_fmri_nii'; 

%% set subject's root directory
n_sbj = 31;


for sbj_i = 1: n_sbj

    % convert subject number start from 1
    n_sbj = strcat("sub-",num2str(sbj_i, '%02.f'));

% raw data folder
raw_sbj_path = fullfile(file_path_in, string(n_sbj), 'study*', 'series_*');
raw_sbj_folders = dir(raw_sbj_path);
for i=1:length(raw_sbj_folders)
fo = fullfile(raw_sbj_folders(i).folder, raw_sbj_folders(i).name);
fi = dir(fullfile(fo, '*.IMA'));  % Find .IMA files
end

%% raw data folder: T1 / field-mag / field-phase / func
raw_T1_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "T1"),1,"last"));
raw_fmap_mag_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"first"));
raw_fmap_phase_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"last"));
raw_func_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "TASK"),1,"last"));

%% make output directory
output_dir = fullfile(file_path_out,n_sbj);
if ~exist(output_dir,"dir")
    mkdir(output_dir);
    mkdir(output_dir+'\anat');
    mkdir(output_dir+'\fmap\mag');
    mkdir(output_dir+'\fmap\phase');
    mkdir(output_dir+'\func');
end
% addpath(genpath(output_dir));C

%% converting T1 image
folder = fullfile(raw_T1_dir(1).folder, raw_T1_dir.name);
files = dir(fullfile(folder, '*.IMA'));  % Find .IMA files
sbj_T1_images = fullfile({files.folder}, {files.name});  % Create full paths for each file.

matlabbatch{1}.spm.util.import.dicom.data = sbj_T1_images(:);
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = cellstr({fullfile(output_dir, 'anat')});
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

%% converting field map - magnitude
folder = fullfile(raw_fmap_mag_dir(1).folder, raw_fmap_mag_dir.name);
files = dir(fullfile(folder, '*.IMA'));  
sbj_fmap_mag = fullfile({files.folder}, {files.name});

matlabbatch{2}.spm.util.import.dicom.data = sbj_fmap_mag(:);
matlabbatch{2}.spm.util.import.dicom.root = 'flat';
matlabbatch{2}.spm.util.import.dicom.outdir = cellstr({fullfile(output_dir, 'fmap', 'mag')});
matlabbatch{2}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{2}.spm.util.import.dicom.convopts.format = 'img';
matlabbatch{2}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{2}.spm.util.import.dicom.convopts.icedims = 0;

%% converting field map - phase
folder = fullfile(raw_fmap_phase_dir(1).folder, raw_fmap_phase_dir.name);
files = dir(fullfile(folder, '*.IMA'));  
sbj_fmap_phase = fullfile({files.folder}, {files.name});

matlabbatch{3}.spm.util.import.dicom.data = sbj_fmap_phase(:);
matlabbatch{3}.spm.util.import.dicom.root = 'flat';
matlabbatch{3}.spm.util.import.dicom.outdir = cellstr({fullfile(output_dir, 'fmap', 'phase')});
matlabbatch{3}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{3}.spm.util.import.dicom.convopts.format = 'img';
matlabbatch{3}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{3}.spm.util.import.dicom.convopts.icedims = 0;

%% converting MR image
folder = fullfile(raw_func_dir(1).folder, raw_func_dir.name);
files = dir(fullfile(folder, '*.IMA'));  
sbj_MR_images = fullfile({files.folder}, {files.name});

matlabbatch{4}.spm.util.import.dicom.data = sbj_MR_images(:);
matlabbatch{4}.spm.util.import.dicom.root = 'flat';
matlabbatch{4}.spm.util.import.dicom.outdir = cellstr({fullfile(output_dir, 'func')});
matlabbatch{4}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{4}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{4}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{4}.spm.util.import.dicom.convopts.icedims = 0;



spm_jobman('run', matlabbatch);

end

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


