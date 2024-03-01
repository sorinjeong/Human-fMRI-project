% 지연이가 알려준대로 1st level에서 개별적으로 뽑은 beta값들을 가져와서 mask==1인 부분만 mean을 낸 후 
% 그 mean이 하나의 dot이 되어서 피험자들 모두를 bar plot 그리기

mni_roi_mask_file = 'D:\fMRI\OCAT_DIR\data\atlas\Standard-space\mni152_hippocampal_area.nii';
mni_hpc_coord=mask2coord(mni_roi_mask_file);



% spm_mat_file = 'D:\fMRI\OCAT_DIR\data\240228\2nd_Level\correct\SPM.mat';
% load(spm_mat_file);
cd('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');

path_root = 'D:\fMRI\OCAT_DIR\data\240228';
path_out_1st = fullfile(path_root,'1st_Level'); 
path_ROI = fullfile(path_root,'ROI');mkdir(path_ROI);
load([path_root '\contrast_v1.mat']);
load('event_regressor_all2');

        
for ct_con=1:numel(contrast{1,2})
    for n_sbj=1:numel(sbj_id_list)
        c_sbj = sprintf('sub-%.2d', sbj_id_list(n_sbj));
        
        sbj_beta = fullfile(path_out_1st, c_sbj, sprintf('%s%04d%s', 'con_',ct_con,'.nii'));
        roi_beta = niftiread(sbj_beta);
        
        A=[];
        for i=1:size(mni_hpc_coord,1)
            x=mni_hpc_coord(i,1); y=mni_hpc_coord(i,2); z=mni_hpc_coord(i,3);
            A(i)=roi_beta(x,y,z);
            roi_mean(n_sbj) = nanmean(A);
        end  
    end
    conditions_roi{ct_con,1} = contrast{1,2}{ct_con};
    conditions_roi{ct_con,2} = roi_mean;
end
        
        
        
        
