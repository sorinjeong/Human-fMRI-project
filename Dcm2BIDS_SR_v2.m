addpath(genpath('C:\spm'));addpath(genpath('D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data'));
% addpath('C:\spm\spm12\external\fieldtrip\compat\matlablt2017b');
%%
raw_data_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\brain_raw';
preprocessed_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\1012practice';

%% set subject's root directory
% sbjNum = [85, 90, 96, 100];
sbjNum = [102];
sbjName = strcat('SUB',string(sbjNum));
raw_sbj_path = fullfile(raw_data_path, sbjName, '**', 'series_*');
raw_sbj_folders = dir(raw_sbj_path);

%% raw data folder: T1 / field-mag / field-phase / func
raw_T1_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "T1")));
raw_fmap_mag_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"first"));
raw_fmap_phase_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "AP"),1,"last"));
raw_func_dir = raw_sbj_folders(find(contains({raw_sbj_folders.name}, "TASK5")));
%% make output directory
output_dir = fullfile(preprocessed_path,sbjName);
mkdir(output_dir);
    mkdir(output_dir+'\anat');
    mkdir(output_dir+'\fmap\mag');
    mkdir(output_dir+'\fmap\phase');
    mkdir(output_dir+'\func');
addpath(genpath(output_dir));

%% converting T1 image
folder = fullfile(raw_T1_dir(1).folder, raw_T1_dir.name);
files = dir(fullfile(folder, '*.IMA'));  % Find .IMA files
sbj_T1_images = fullfile({files.folder}, {files.name});  % Create full paths for each file.

matlabbatch{1}.spm.util.import.dicom.data = sbj_T1_images(:);
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = {(output_dir+'\anat')};
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
matlabbatch{2}.spm.util.import.dicom.outdir = {(output_dir+'\fmap\mag')};
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
matlabbatch{3}.spm.util.import.dicom.outdir = {(output_dir+'\fmap\phase')};
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
matlabbatch{4}.spm.util.import.dicom.outdir = {(output_dir+'\func')};
matlabbatch{4}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{4}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{4}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{4}.spm.util.import.dicom.convopts.icedims = 0;


spm('defaults','FMRI');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);


