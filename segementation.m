addpath('C:\Users\Leelab\Documents\MATLAB\spm12');
path_in='D:\MJ\data_fmri_organized_seg_final\data_fmri_organized_seg_final';
cd(path_in)
load('regressors_GLM_0717.mat')
sbj_id_list_38(sbj_id_list_38==7)=[];

% bi_HP_variable='bi_HP';

for sbj_i=1:34
    n_sbj=sbj_id_list_38(sbj_i);
    c_sbj=sprintf('sub-%.2d',n_sbj);

path_out=strcat('./nii/',c_sbj);

mkdir(path_out)


load(fullfile(path_in,strcat(string(n_sbj),'.mat')))


header=seg.seg_fit.header;
header.Datatype='double';
for i=1:length(seg.seg_fit.hpc.roi_list_nn)
roi_mask=double(seg.seg_fit.hpc.roi_list_nn{i});
roi_name=seg.seg_fit.hpc.roi_name_list{i};
% output_filename=sprintf('sub-%.2d./nii/',bi_HP_variable,i);
output_filename=sprintf('%s.nii',roi_name);

niftiwrite(roi_mask,fullfile(path_out,output_filename),header);
end
end

