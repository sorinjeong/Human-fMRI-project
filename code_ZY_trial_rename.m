%% single trial naming
zcode_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
addpath(zcode_path)
load("sbj_events.mat")
singletrial_path='C:\Users\User\Desktop\JSR\GLMlocal_v5\1stLevelSinggleTrial';
addpath(singletrial_path)
load("sbj_id_list.mat")
sbj_id_list(sbj_id_list==7)=[];




for sbj_i=1:numel(sbj_id_list)
    sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));

    sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001','OC_*'));

    for i=1:height(sbj_OC)

        if contains(sbj_OC(i).name,'OC_asso_obj_corr')
            roc='hit';
        elseif contains(sbj_OC(i).name,'OC_asso_obj_incorr')
            roc='miss';
        elseif contains(sbj_OC(i).name,'OC_NA_obj_corr')
            roc='corr_rej';
        elseif contains(sbj_OC(i).name,'OC_NA_obj_incorr')
            roc='false';
        end

        num=str2double(sbj_OC(i).name(end-1:end));
        match_cond=find(tbl.session == sbj_id_list(sbj_i) & strcmp(tbl.roc, roc));
        trial_n=tbl.Trial(match_cond(num));


beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
beta_rename=strcat(string(trial_n),'_',roc,'_',sbj_OC(i).name,'.nii');

copyfile(beta, fullfile(singletrial_path,sbj_n, beta_rename));


    end
end