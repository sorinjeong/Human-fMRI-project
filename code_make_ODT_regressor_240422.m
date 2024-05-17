
clc;clear;

root_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';cd(root_path)
load('sbj_id_list')
%%
load('regressors_GLM_0513');

reg_for_glm=reg_for_glm_ver1;

for i = 1:length(sbj_id_list)
    O_names = reshape(reg_for_glm{1,i}.ODT_detail.regress_name',[],1);
    O_onsets = reshape(reg_for_glm{1,i}.ODT_detail.regress_onset',[],1);
    O_durations = reshape(reg_for_glm{1,i}.ODT_detail.regress_duration',[],1);

    %% create reg_for_glm
    odt_reg{1, i}.ODT= reg_for_glm_ver1{1,i}.ODT_detail;
    odt_reg{1, i}.ODT.regress_name = O_names(1:10);
    odt_reg{1, i}.ODT.regress_onset = O_onsets(1:10);
    odt_reg{1, i}.ODT.regress_duration = O_durations(1:10);

    reg_for_glm_ver1{1,i}.ODT_obj4=odt_reg{1, i}.ODT;
    reg_for_glm_ver2{1,i}.ODT_obj4=odt_reg{1, i}.ODT;

end
save(fullfile('regressors_ODT_4obj_0513'),"odt_reg","sbj_id_list");

save(fullfile('regressors_GLM_0513'),"odt_reg","reg_for_glm_ver1","reg_for_glm_ver2","sbj_id_list");





