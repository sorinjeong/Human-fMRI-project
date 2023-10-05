raw_data_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\brain_raw';
preprocessed_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\1004practice\1st';

%% set subject's root directory
% sbjNum = [85, 90, 96, 100];
sbjNum = [102];
sbjName = strcat('SUB',string(sbjNum));
raw_sbj_path = fullfile(raw_data_path, sbjName, 'study*', 'series_*');
raw_sbj_folder = dir(raw_sbj_path);
% raw_sbj_dir = fullfile(raw_sbj_path, {raw_sbj_folder.name});


%% raw data folder: T1 / field-mag / field-phase / func
raw_T1_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "T1")));
raw_fmap_mag_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "AP"),1,"first"));
raw_fmap_phase_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "AP"),1,"last"));
raw_func_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "TASK5")));



%% make output directory
output_dir = fullfile(preprocessed_path, sbjName);
if ~exist("output_dir")
    mkdir(fullfile(output_dir, "anat"))
    mkdir(fullfile(output_dir, "fmap","mag"))
    mkdir(fullfile(output_dir, "fmap","phase"))
    mkdir(fullfile(output_dir, "func"))
end

%% converting T1 image
sbj_T1_images = {dir(fullfile(raw_T1_folder(1).folder, raw_T1_folder.name)).name}';

matlabbatch{1}.spm.util.import.dicom.data = sbj_T1_images;
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = {fullfile(output_dir, "anat")};
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

%% converting field map - magnitude
sbj_fmap_mag = {dir(fullfile(raw_fmap_mag_folder(1).folder, raw_fmap_mag_folder.name)).name}';

matlabbatch{2}.spm.util.import.dicom.data = sbj_fmap_mag;
matlabbatch{2}.spm.util.import.dicom.root = 'flat';
matlabbatch{2}.spm.util.import.dicom.outdir = {fullfile(output_dir, "fmap","mag")};
matlabbatch{2}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{2}.spm.util.import.dicom.convopts.format = 'img';
matlabbatch{2}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{2}.spm.util.import.dicom.convopts.icedims = 0;

%% converting field map - phase
sbj_fmap_phase = {dir(fullfile(raw_fmap_phase_folder(1).folder, raw_fmap_phase_folder.name)).name}';

matlabbatch{3}.spm.util.import.dicom.data = sbj_fmap_phase;
matlabbatch{3}.spm.util.import.dicom.root = 'flat';
matlabbatch{3}.spm.util.import.dicom.outdir = {fullfile(output_dir, "fmap","phase")};
matlabbatch{3}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{3}.spm.util.import.dicom.convopts.format = 'img';
matlabbatch{3}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{3}.spm.util.import.dicom.convopts.icedims = 0;

%% converting MR image
sbj_MR_images = {dir(fullfile(raw_func_folder(1).folder, raw_func_folder.name)).name}';

matlabbatch{4}.spm.util.import.dicom.data = sbj_MR_images;
matlabbatch{4}.spm.util.import.dicom.root = 'flat';
matlabbatch{4}.spm.util.import.dicom.outdir = {fullfile(output_dir, "func")};
matlabbatch{4}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{4}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{4}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{4}.spm.util.import.dicom.convopts.icedims = 0;


    batch = matlabbatch;
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);
    clear batch