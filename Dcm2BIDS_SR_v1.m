raw_data_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\brain_raw';
preprocessed_path = 'D:\leelab\Human fMRI projects\0O-CAT\fMRI_ocat_data\1004practice\1st';

%% set subject's root directory
% sbjNum = [85, 90, 96, 100];
sbjNum = [102];
sbjName = strcat('SUB',string(sbjNum));
raw_sbj_path = fullfile(raw_data_path, sbjName, '**', 'series_*');
raw_sbj_folder = dir(raw_sbj_path);
% raw_sbj_dir = fullfile(raw_sbj_path, {raw_sbj_folder.name});


%% raw data folder: T1 / field-mag / field-phase / func
raw_T1_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "T1")));
raw_fmap_mag_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "AP"),1,"first"));
raw_fmap_phase_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "AP"),1,"last"));
raw_func_folder = raw_sbj_folder(find(contains({raw_sbj_folder.name}, "TASK5")));



%% make output directory
output_dir = [preprocessed_path sbjName];
if exist("output_dir")
    mkdir([output_dir '\anat'])
    mkdir([output_dir '\fmap\mag'])
    mkdir([output_dir '\fmap\phase'])
    mkdir([output_dir '\func'])
end


%% converting T1 image
sbj_T1_images = dir(fullfile(raw_T1_folder(1).folder, raw_T1_folder.name)).name;
matlabbatch{1}.spm.util.import.dicom.data = sbj_T1_images;





