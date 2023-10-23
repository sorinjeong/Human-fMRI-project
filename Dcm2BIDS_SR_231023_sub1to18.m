clear all; clc;
restoredefaultpath; % Restore the MATLAB path to its default state.
addpath('C:\Users\Leelab\Documents\MATLAB\spm12');
addpath('C:\Users\Leelab\Documents\MATLAB\spm12\external\fieldtrip'); % Add the FieldTrip directory back to the MATLAB path.
ft_defaults; % Reset the FieldTrip defaults.

addpath(genpath('\\147.47.203.148\LeeStorage2\EPhysRawData\fmri_ocat_analysis\O-CAT_fMRIdata_backup'));

%%
spm('defaults','FMRI');
spm_jobman('initcfg');

%%
raw_data_path = '\\147.47.203.148\LeeStorage2\EPhysRawData\fmri_ocat_analysis\O-CAT_fMRIdata_backup\brain_raw';
preprocessed_path = 'C:\Users\leelab\Desktop\fMRI\OCAT\OCAT_DIR';

%% set subject's root directory
sbjNum = [85:102];
% sbjNum = [102];
sbjName = strcat('SUB',string(sbjNum));

for sub = 1: length(sbjName)
    % convert subject number start from 1
    if sub<10
    sbj = strcat("sub-0",string(sub)); 
    else
    sbj = strcat("sub-",string(sub));
    end
% raw data folder
raw_sbj_path = fullfile(raw_data_path, sbjName(sub), 'study*', 'series_*');
raw_sbj_folders = dir(raw_sbj_path);

%% raw data folder: T1 / field-mag / field-phase / func
raw_T1_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "T1")));
raw_fmap_mag_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"first"));
raw_fmap_phase_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"last"));
raw_func_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "TASK5")));

%% make output directory
output_dir = fullfile(preprocessed_path,sbj);
if ~exist(output_dir,"dir")
    mkdir(output_dir);
    mkdir(output_dir+'\anat');
    mkdir(output_dir+'\fmap\mag');
    mkdir(output_dir+'\fmap\phase');
    mkdir(output_dir+'\func');
end
% addpath(genpath(output_dir));

%% converting T1 image
folder = fullfile(raw_T1_dir(1).folder, raw_T1_dir.name);
files = dir(fullfile(folder, '*.IMA'));  % Find .IMA files
if ~isempty(dir(fullfile(folder, 'SUB8*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB9*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB10*.IMA')))
for idx = 1:numel(files)
    movefile(fullfile(folder, files(idx).name), fullfile(folder, strrep(files(idx).name, string(sbjName(sub)), string(sbj))));
end
end
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
if ~isempty(dir(fullfile(folder, 'SUB8*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB9*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB10*.IMA'))) 
for idx = 1:numel(files)
    movefile(fullfile(folder, files(idx).name), fullfile(folder, strrep(files(idx).name, string(sbjName(sub)), string(sbj))));
end
end
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
if ~isempty(dir(fullfile(folder, 'SUB8*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB9*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB10*.IMA')))
for idx = 1:numel(files)
    movefile(fullfile(folder, files(idx).name), fullfile(folder, strrep(files(idx).name, string(sbjName(sub)), string(sbj))));
end
end
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
if ~isempty(dir(fullfile(folder, 'SUB8*.IMA'))) || ~isempty(dir(fullfile(folder, 'SUB9*.IMA')))|| ~isempty(dir(fullfile(folder, 'SUB10*.IMA')))
for idx = 1:numel(files)
    movefile(fullfile(folder, files(idx).name), fullfile(folder, strrep(files(idx).name, string(sbjName(sub)), string(sbj))));
end
end
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
