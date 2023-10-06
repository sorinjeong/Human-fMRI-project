raw_data_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\brain_raw';
preprocessed_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\1004practice\1st';

%% set subject's root directory
% sbjNum = [85, 90, 96, 100];
sbjNum = [102];
sbjName = strcat('SUB',string(sbjNum));
raw_sbj_path = dir(fullfile(raw_data_path, sbjName, 'study*', 'series_*'));
raw_sbj_folder = {raw_sbj_path.name};
raw_sbj_path = raw_sbj_path.folder;
% raw_sbj_dir = fullfile(raw_sbj_path, {raw_sbj_folder.name});


%% raw data folder: T1 / field-mag / field-phase / func
raw_T1_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "T1")));
raw_fmap_mag_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "AP"),1,"first"));
raw_fmap_phase_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "AP"),1,"last"));
raw_func_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "TASK5")));

dir(fullfile(raw_sbj_path.folder,'*T1*'))


category = ["T1" "AP" "TASK5"];rawdata_folders=[];
for i=1:length(category)
foldername = raw_sbj_folder(find(contains(raw_sbj_folder,category(i))));
if i==1
raw_T1_folder = string(strcat(raw_sbj_path, '\', foldername));
raw_T1_files = dir(raw_T1_folder).name
raw_T1_filelist = 

dir(['*.IMA'])

d= uigetdir(raw_sbj_path, foldername)
uigetfile('*.IMA', 'Select Multiple Files', raw_T1_folder, 'MultiSelect','on')

if i==2
rawdata_folders= [rawdata_folders; foldername(1); foldername(2)];
else
rawdata_folders= [rawdata_folders; foldername];
end
end


%% make output directory
output_dir = fullfile(preprocessed_path, sbjName);
if ~exist("output_dir")
    mkdir(fullfile(output_dir, "anat"))
    mkdir(fullfile(output_dir, "fmap","mag"))
    mkdir(fullfile(output_dir, "fmap","phase"))
    mkdir(fullfile(output_dir, "func"))
end

%% converting T1 image
sbj_T1_images = {dir(fullfile(raw_T1_folder.folder, raw_T1_folder.name)).name};

matlabbatch{1}.spm.util.import.dicom.data = sbj_T1_images;
matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = {fullfile(output_dir, "anat")};
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

%% converting field map - magnitude
sbj_fmap_mag = {dir(fullfile(raw_fmap_mag_folder.folder, raw_fmap_mag_folder.name)).name};

matlabbatch{2}.spm.util.import.dicom.data = sbj_fmap_mag;
matlabbatch{2}.spm.util.import.dicom.root = 'flat';
matlabbatch{2}.spm.util.import.dicom.outdir = {fullfile(output_dir, "fmap","mag")};
matlabbatch{2}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{2}.spm.util.import.dicom.convopts.format = 'img';
matlabbatch{2}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{2}.spm.util.import.dicom.convopts.icedims = 0;

%% converting field map - phase
sbj_fmap_phase = {dir(fullfile(raw_fmap_phase_folder.folder, raw_fmap_phase_folder.name)).name};

matlabbatch{3}.spm.util.import.dicom.data = sbj_fmap_phase;
matlabbatch{3}.spm.util.import.dicom.root = 'flat';
matlabbatch{3}.spm.util.import.dicom.outdir = {fullfile(output_dir, "fmap","phase")};
matlabbatch{3}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{3}.spm.util.import.dicom.convopts.format = 'img';
matlabbatch{3}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{3}.spm.util.import.dicom.convopts.icedims = 0;

%% converting MR image
sbj_MR_images = {dir(fullfile(raw_func_folder.folder, raw_func_folder.name)).name};

matlabbatch{4}.spm.util.import.dicom.data = sbj_MR_images;
matlabbatch{4}.spm.util.import.dicom.root = 'flat';
matlabbatch{4}.spm.util.import.dicom.outdir = {fullfile(output_dir, "func")};
matlabbatch{4}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{4}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{4}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{4}.spm.util.import.dicom.convopts.icedims = 0;


    batch = matlabbatch;
    spm('defaults','FMRI');
    spm_jobman('initcfg');
    spm_jobman('run',batch);
    clear batch



tbl = struct2table(raw_sbj_folder)
tbl={horzcat(tbl.folder,tbl.name)}
tbl(:,[2 1])
mergevars(tbl,[2 1])



