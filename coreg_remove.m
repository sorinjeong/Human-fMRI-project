clear;clc;

cd 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code'
addpath('C:\Users\User\Documents\MATLAB\');
load('regressors_GLM_0405');
% sub-10 제외.
sbj_id_list(sbj_id_list==10)=[];
%sub-07 추가
sbj_id_list=sort([sbj_id_list 7]);

TYPE='HBT'; % 'HBT' or 'FS60'
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\' TYPE];


roi=struct;
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'coreg_coreg_*.nii'));
   for i = 1:length(local_masks)
    old_name = local_masks(i).name;
    new_name = strrep(old_name, 'coreg_', '');
    movefile(fullfile(mask_path, sprintf('s-%.2d',sbj_n), old_name), fullfile(mask_path, sprintf('s-%.2d',sbj_n), new_name));
end
end