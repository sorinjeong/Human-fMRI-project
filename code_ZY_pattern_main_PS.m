%% single trial naming
zcode_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
addpath(zcode_path)
load("sbj_events.mat")
singletrial_path='C:\Users\User\Desktop\JSR\GLMlocal_v5\1stLevelSinggleTrial'; %90번
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

if trial<10
beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
beta_rename=strcat(sprintf('%.2d',trial_n),'_',roc,'_',sbj_OC(i).name,'.nii');

copyfile(beta, fullfile(singletrial_path,sbj_n, beta_rename));

end
    end
end



%% Trial pattern extraction
% % trial별 roc name (hit, miss, corr_rej, false)
% load("sbj_events.mat")
% 
% roc_cell=cell(1, ceil(length(roc_trial)/32));
% for i = 1:length(roc_cell)
%     startIdx = (i-1)*32 + 1;
%     endIdx = min(i*32, length(roc_trial));
%     roc_cell{i} = roc_trial(startIdx:endIdx);
% end

% % trial별 object 및 context 정보
% obj_num=cell(1, ceil(height(tbl)/32));
% context=cell(1, ceil(height(tbl)/32));
% for i = 1:length(obj_num)
%     startIdx = (i-1)*32 + 1;
%     endIdx = min(i*32, height(tbl));
%     obj_num{i} = tbl.Object(startIdx:endIdx)';
%     for cont=1:32
%         if tbl.Association(cont) ==1
%             ctxt{cont} = tbl.Context(cont);
%         elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==1
%             ctxt{cont} = 2;
%         elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==2
%             ctxt{cont} = 1;
%         end
%         context{i}=ctxt;
%     end
% end


% pattern structure
roi_mat_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\atlas\ZY_240612\data_fmri_organized_seg_v5_SR';

OC_pattern=struct;
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    load(fullfile(roi_mat_path,string(sbj_n)))
    roi=seg.seg_fit.hpc;
    sbj_nii=dir(fullfile(singletrial_path,sprintf('sub_%d',sbj_n),'*.nii'));
    curr_roi=roi.roi_list_nn;
    for r=1:numel(roi.roi_name_list)
        for t=1:32
            curr_beta=niftiread(fullfile(sbj_nii(1).folder,sbj_nii(t).name));

            pattern_roi=curr_beta(curr_roi{1,r});
            OC_pattern.data{1,sbj_i}{1,r}{1,t}=pattern_roi;
            OC_pattern.data{1,sbj_i}{1,r}{2,t}=roc_cell{1,sbj_i}{t};
            OC_pattern.data{1,sbj_i}{1,r}{3,t}=obj_num{1,sbj_i}(t);%object정보
            OC_pattern.data{1,sbj_i}{1,r}{4,t}=context{1,sbj_i}{t}; %context정보
        end

    end
    OC_pattern.roi=roi;

end
% subject(27) - roi region(30) - trial(32)


save('C:\Users\User\Desktop\JSR\GLMlocal_v5\OC_pattern.mat',"OC_pattern");





%% same Ps, diff PS
parsing=struct;
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    for r=1:numel(roi.roi_name_list)

        curr_s_r=OC_pattern.data{1,sbj_i}{1,r};
        bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');


        f_idx = find(cellfun(@(x) x == 1, curr_s_r(4,1:4)));
        curr_f= cell2mat(OC_pattern.data{1,sbj_i}{1,r}(3,f_idx));
        curr_c= setdiff(4:7,curr_f);

        f1_idx=find(cell2mat(curr_s_r(3,:))==curr_f(1));
        f2_idx=find(cell2mat(curr_s_r(3,:))==curr_f(2));
        c1_idx=find(cell2mat(curr_s_r(3,:))==curr_c(1));
        c2_idx=find(cell2mat(curr_s_r(3,:))==curr_c(2));

%% same
same_ctxt.all{1,sbj_i}=[bulk_corr(f1_idx,f2_idx); bulk_corr(c1_idx,c2_idx)];

m=[];
for lap=1:8
    m(lap)=mean([same_ctxt.all{1,sbj_i}(lap,lap); same_ctxt.all{1,sbj_i}(lap+8,lap)],"all");
end
same_ctxt.lap{1,sbj_i}{:,r}=m';
for block=1:4
    same_ctxt.block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
end

%% diff
% complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2

diff_ctxt.all{1,sbj_i}{1}=bulk_corr(f1_idx,c1_idx);
diff_ctxt.all{1,sbj_i}{2}=bulk_corr(f1_idx,c2_idx);
diff_ctxt.all{1,sbj_i}{3}=bulk_corr(f2_idx,c1_idx);
diff_ctxt.all{1,sbj_i}{4}=bulk_corr(f2_idx,c2_idx);

m=[];
for lap=1:8
    m(lap)=mean([diff_ctxt.all{1,sbj_i}{1}(lap,lap); diff_ctxt.all{1,sbj_i}{2}(lap,lap);diff_ctxt.all{1,sbj_i}{3}(lap,lap);diff_ctxt.all{1,sbj_i}{4}(lap,lap)],"all");
end
diff_ctxt.lap{1,sbj_i}{:,r}=m';
for block=1:4
    diff_ctxt.block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
end



    end
    same_ctxt.half{1,sbj_i}(1,:)=mean(cell2mat(same_ctxt.block{1,sbj_i}(1:2,:)));
    same_ctxt.half{1,sbj_i}(2,:)=mean(cell2mat(same_ctxt.block{1,sbj_i}(3:4,:)));
    diff_ctxt.half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.block{1,sbj_i}(1:2,:)));
    diff_ctxt.half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.block{1,sbj_i}(3:4,:)));
end







%% 지연 PS 랑 맞추기!!!

fail_sbj=[5 6 15 22 30];

same=struct;diff=struct;

which_sbj='learn'; %'all', 'learn', 'fail'
roi = {'Lt.Hp', 'Lt.PHC', 'Lt.CA23DG', 'Lt.CA1', ...
    'Rt.Hp', 'Rt.PHC',  'Rt.CA23DG','Rt.CA1', ...
    'Bi.Hp', 'Bi.PHC',  'Bi.CA23DG','Bi.CA1','Bi.EC'};
% region indexing
region_idx=find(ismember(hpc_name_all, roi));
hpc_region=hpc_name_all(region_idx);

% sbj indexing
if strcmp(which_sbj,'learn')
    sbj_idx= find(~ismember(sbj_id_list, fail_sbj));
    ttl = "learned subjects (n=22)";

elseif strcmp(which_sbj,'fail')
    sbj_idx= find(ismember(sbj_id_list, fail_sbj));
    ttl = "failed subjects (n=5)";
else
    sbj_idx= 1:numel(sbj_id_list);
    ttl = "All subjects(n=27)";
end



for sbj_i=1:numel(sbj_idx)
    sbj_n=sbj_idx(sbj_i);
    %     for bl=1:4
    for rg=1:numel(hpc_region)
        rg_n=strrep(hpc_region{rg},'.','_');
        %         same.(which_sbj).(sprintf('block%d',bl))(sbj_i,:)=same_category_all{1, sbj_n}(bl,region_idx);
        %         diff.(which_sbj).(sprintf('block%d',bl))(sbj_i,:)=diff_asso_na_all{1, sbj_n}(bl,region_idx);
        sameZY.(which_sbj).(rg_n)(sbj_i,:)=same_category_all{1, sbj_n}(:,region_idx(rg))';
        diffZY.(which_sbj).(rg_n)(sbj_i,:)=diff_asso_na_all{1, sbj_n}(:,region_idx(rg))';
        sameSR.(which_sbj).(rg_n)(sbj_i,:)=same_ctxt.block{1, sbj_n}(:,region_idx(rg))';
        diffSR.(which_sbj).(rg_n)(sbj_i,:)=diff_ctxt.block{1, sbj_n}(:,region_idx(rg))';


        same_half_ZY.(rg_n)(:,1)=mean(sameZY.(which_sbj).(rg_n)(:,1:2),2);
        same_half_ZY.(rg_n)(:,2)=mean(sameZY.(which_sbj).(rg_n)(:,3:4),2);

        diff_half_ZY.(rg_n)(:,1)=mean(diffZY.(which_sbj).(rg_n)(:,1:2),2);
        diff_half_ZY.(rg_n)(:,2)=mean(diffZY.(which_sbj).(rg_n)(:,3:4),2);
    end
end


%--> 이거로 69번 컴에서 graph pad prism 으로 plot그린거임 (block별, half 별)











































