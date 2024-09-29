% %% single trial naming
% % clc;clear
% 
% % zcode_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
% % addpath(zcode_path)
% % root_path = 'C:\Users\User\Desktop\JSR';%90번컴
% root_path = 'G:\JSR\240922_new_fmri'; %69번컴
% addpath(fullfile(root_path, 'spm_prep_glm_0922')); addpath('G:\JSR\sbj39to41')
% cd(fullfile(root_path, 'spm_prep_glm_0922','main'))
% % zroot_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\spm_prep_glm_0725';
% % addpath(zroot_path)
% 
% load('../../bhv/data_bhv_log_table/total/num_sbj_events.mat')
% load('regressors_GLM_0922_backup.mat')
% load('../../bhv/learn_curv_info.mat');
% roi_hpc_name=strrep(roi_hpc_name,'CA23DG','DGCA3');
% 
% % sbj_id_list=setdiff(sbj_id_list_41,[2,3,5,6,7,15,21,22,30,32,33,38]);
% sbj_id_list=sbj_id_list_41;
% tbl=num_sbj_events;
% % setting
% main_or_ODT = {'main'}; % {'main','pre-ODT','post-ODT'} for for loops
% 
% 
% 
% hpc_select=[1:5,7:17,19:24];
% ctx_select=1:30;
% 
% %%
% 
% for m_o = main_or_ODT
%     disp(m_o)
% 
%         singletrial_path = fullfile(root_path, 'spm_prep_glm_0922', m_o{:}, 'single_trial');
%     
%     %     for sbj_i=1:numel(sbj_id_list)
%     %                 sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
%     %     files=dir(fullfile(singletrial_path,sbj_n, '*_0*'));
%     %     for i=1:height(files)
%     %     delete(fullfile(files(1).folder,files(i).name))
%     %     end
%     %     end
%     %     end
%     
%     %% rename_ trial
%     if strcmp(m_o{:}, 'main')
%         disp('rename_ trial')
%         
%         %         singletrial_path = fullfile(zroot_path, m_o{:},ver, 'single_trial');
%         
%         for sbj_i=1:numel(sbj_id_list)
%             sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
%             %             singletrial_path_out=fullfile(root_path, 'spm_prep_glm_0922', m_o{:},ver, 'single_trial',sbj_n);
%             %             mkdir(singletrial_path_out)
%             
%             
%             pp='obj_show';
%             sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001',strcat(pp,'_*')));
%             roc=cellfun(@(x) extractBetween(x,strcat(pp,'_'),'_0'),{sbj_OC.name});
%             d=reg_0922{1, sbj_i}.trial_detail;
%             un=unique(roc);
%             
%             sidx=[];
%             for r=1:numel(un)
%                 idx=d.(un{r})';
%                 sidx=[sidx idx];
%             end
%             [b,iii]=sort(sidx);
%             if max(iii)==32
%                 tbl.roc(tbl.session==sbj_id_list(sbj_i))=roc(iii);
%                 for i=1:height(sbj_OC)
%                     beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
%                     beta_rename=sprintf('%.2d_%s.nii',sidx(i), sbj_OC(i).name);
%                     copyfile(beta, fullfile(singletrial_path,sbj_n, beta_rename));
%                 end
%             else
%                 error('The total number of trials is not 32!');
%             end
%         end
%         
%     else
%         %         singletrial_path = fullfile(root_path, 'spm_prep_glm_0922', m_o{:}, 'single_trial');
%         for sbj_i=1:numel(sbj_id_list)
%             sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
%             pp=extractBefore(m_o{:},'-');
%             sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001',strcat(pp,'*')));
%             obj=cellfun(@(x) extractBetween(x,strcat(pp,'_'),'_0'),{sbj_OC.name},'UniformOutput',false);
%             %         obj=repelem(["f1","f2","c1","c2","target"],3);
%             temp_onset=reg_0922{1, sbj_i}.ODT.(pp).regress_onset;
%             if contains(obj{1},'city')
%                 temp_onset=[temp_onset(3:4); temp_onset(1:2); temp_onset(5)]; %["f1","f2","c1","c2","target"] regressor 순서로 바꿔줌
%             end
%             [tr,tr_idx]=sort(cell2mat(temp_onset));
%             
%             %    temp=table2array(readtable('Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\0716\data_bhv_log_table_38\total\event_table_MR_new.xlsx','Sheet',sprintf('sub-%d',sbj_id_list(sbj_i)),'Range', 'F:F'));
%             %
%             %       ODT_context.pre{sbj_i}=temp(1:15);
%             %     ODT_context.post{sbj_i}=temp(48:end);
%             %
%             
%             %             beta_rename={};
%             %             for tridx=1:length(tr_idx)
%             %                 beta_rename{tridx}=sprintf('%.2d_%s',tridx,sbj_OC(tr_idx(tridx)).name);
%             %             end
%             %             for i=1:numel(beta_rename)
%             %                 beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
%             %                 copyfile(beta, fullfile(singletrial_path,sbj_n, strcat(beta_rename{i},'.nii')));
%             %             end
%             %
%             
%             
%             beta_rename={};
%             for tridx=1:length(tr_idx)
%                 beta_rename{tridx}=sprintf('%.2d_%s',tridx,sbj_OC(tr_idx(tridx)).name);
%                 beta=fullfile(sbj_OC(tr_idx(tridx)).folder,sbj_OC(tr_idx(tridx)).name,'beta_0001.nii');
%                 copyfile(beta, fullfile(singletrial_path,sbj_n, strcat(beta_rename{tridx},'.nii')));
%             end
%             
%             
%             reg_0922{1, sbj_i}.ODT.(pp).object=obj(tr_idx);
%             reg_0922{1, sbj_i}.ODT.(pp).object_idx=tr_idx;
%         end
%         save('regressors_GLM_0922.mat',"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name')
%     end
%     %     save(string(fullfile(root_path, 'spm_prep_glm_0922','regressors_GLM_0922.mat')),"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name')
%     
%     
%     
%     %% Trial detail
%     % load("sbj_events.mat")
%     
%     if strcmp(m_o{:}, 'main')
%         disp('trial detail')
%         
%         roc_cell=cell(1, ceil(length(tbl.roc)/32));
%         for i = 1:length(roc_cell)
%             startIdx = (i-1)*32 + 1;
%             endIdx = min(i*32, length(tbl.roc));
%             roc_cell{i} = tbl.roc(startIdx:endIdx);
%         end
%         
%         % trial별 object 및 context 정보
%         obj_num=cell(1, numel(sbj_id_list));
%         context=cell(1, numel(sbj_id_list));
%         for sbj_i=1:numel(sbj_id_list)
%             
%             startIdx = (sbj_id_list(sbj_i)-1)*32 + 1;
%             endIdx = min(sbj_id_list(sbj_i)*32, height(tbl));
%             obj_num{sbj_i} = tbl.Object(startIdx:endIdx)';
%             for cont=1:32
%                 if tbl.Association(cont) ==1
%                     ctxt{cont} = tbl.Context(cont);
%                 elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==1
%                     ctxt{cont} = 2;
%                 elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==2
%                     ctxt{cont} = 1;
%                 end
%                 context{sbj_i}=ctxt;
%             end
%             
%             reg_0922{1, sbj_i}.trial_detail.context=cell2mat(context{sbj_i});
%             reg_0922{1, sbj_i}.trial_detail.object=obj_num{sbj_i};
%         end
%         save('regressors_GLM_0922.mat',"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name')
%         
%     end
%     
%     
%     %% pattern structure
%     disp('pattern structure')
%     % load('roi_name.mat')
    
    OC_pattern=struct;

        singletrial_path = fullfile(root_path, 'spm_prep_glm_0922', m_o{:}, 'single_trial');
   
    
    for sbj_i=1:numel(sbj_id_list)
        sbj_n = sbj_id_list(sbj_i);disp(sbj_n)
        load(fullfile(root_path,'data_fmri_seg',string(strcat(string(sbj_n), '.mat'))))
        %         if sbj_i==1
        %             roi_hpc_name=seg.seg_fit.hpc.roi_name_list;
        %             roi_ctx_name=seg.seg_fit.ctx.roi_name_list;
        %         end
        % set roi to number
        for hpc_or_ctx={'hpc','ctx'}
            roi_nn=seg.seg_fit.(hpc_or_ctx{:}).roi_list_nn;
            
            % set roi info
            if strcmp(hpc_or_ctx,'hpc')
                curr_mask=roi_nn(hpc_select);
                curr_names=roi_hpc_name(hpc_select);
                
            else
                curr_mask=roi_nn(ctx_select);
                curr_names=roi_ctx_name(ctx_select);
            end
            sbj_nii=dir(fullfile(singletrial_path,sprintf('sub_%d',sbj_n),'*.nii'));
            
            for r=1:numel(curr_names)
                for t=1:height(sbj_nii) % 32 trials for main / 15 trials for ODT
                    curr_beta=niftiread(fullfile(sbj_nii(1).folder,sbj_nii(t).name));
                    pattern_roi=curr_beta(curr_mask{r});
                    
                    OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{1,t}=pattern_roi;
                    OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{2,t}=extractBetween(sbj_nii(t).name,strcat(pp,'_'),'_0');
                    if strcmp(m_o{:}, 'main')
                        OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,t}=obj_num{1,sbj_i}(t);%object정보
                        OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{4,t}=context{1,sbj_i}{t}; %context정보
                    else
                        OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,t}=reg_0922{1, sbj_i}.ODT.(pp).object{t};
                    end
                end
                % if ~strcmp(m_o{:}, 'main'); OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,:}=reg_0922{1, sbj_i}.ODT.(pp).object(1:12)'; end
                OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{2,r}=curr_names{r};
                
            end
            OC_pattern.(hpc_or_ctx{:}).raw_roi=roi_nn;
        end
        % subject(27) - roi region(30) - trial(32)
    end
    
        save(string(fullfile(root_path,'spm_prep_glm_0922',m_o,'OC_pattern.mat')),"OC_pattern");
    
end
% roi_ctx_name=strrep(roi_ctx_name,'.','_');
% roi_hpc_name=strrep(roi_hpc_name,'.','_');
% save('roi_name.mat','roi_ctx_name','roi_hpc_name')
% 
% %% rewarding
% 
% % %%%%%%%%%%%%%%%%%%% 이 아래 OC pattern에 맞춰서 다시 수정하기!!!
% %% rewarding!!!!
% main_or_ODT = {'main'};
% load('roi_name.mat')
% load('regressors_GLM_0804.mat')

m_o = main_or_ODT;
disp(m_o)
% if strcmp(m_o{:},'main')
rewarding=struct;

% put all/half/block
%     temporal='half1';

for hpc_or_ctx={'hpc','ctx'}
    load(string(fullfile(root_path, 'spm_prep_glm_0922',m_o,'OC_pattern.mat')));
    
    if strcmp(hpc_or_ctx,'hpc')
        curr_names=roi_hpc_name(hpc_select);
    else
        curr_names=roi_ctx_name(ctx_select);
    end
    bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
    
    
    for sbj_i=1:numel(sbj_id_list)
        sbj_n = sbj_id_list(sbj_i);
        for r=1:numel(curr_names)
            
            curr_s_r=OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r};
            bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');
            bulk_corr_ori=bulk_corr;
            
            corr_rej_idx=find(cellfun(@(x) strcmp(x,'corr_rej'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
            hit_idx=find(cellfun(@(x) strcmp(x,'hit'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
            false_idx=find(cellfun(@(x) strcmp(x,'false'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
            miss_idx=find(cellfun(@(x) strcmp(x,'miss'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
            
            %% 진짜 corr_rej_idx만 할건지(correct only), false positive도 섞어서 할건지(corr+incorr).
            % hit+miss / correct rejection + false positive (corr+incorr)로 할거면 아래 코드 실행하기
            corr_rej_idx=[corr_rej_idx false_idx];
            hit_idx=[hit_idx miss_idx];
            
            
            %% all? / half? / block?
            
            % 1. all trials
            % put nothing
            
            % 2. half 1,2
            % half1
                            corr_rej_idx = corr_rej_idx(corr_rej_idx<17);
                            hit_idx = hit_idx(hit_idx<17);
            
            %half2
%             corr_rej_idx = corr_rej_idx(corr_rej_idx>16);
%             hit_idx = hit_idx(hit_idx>16);
            
            
            % 3. block 1,2,3,4
            % block 1
            
            % block 2
            
            % block 3
            
            % block 4
            
            
            
            
            
            rewarding.hit.(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(hit_idx,hit_idx),1);
            rewarding.corr_rej.(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(corr_rej_idx,corr_rej_idx),1);
            
            rewarding.miss.(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(miss_idx,miss_idx),1);
            rewarding.false.(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(false_idx,false_idx),1);
            
            
            
            
            
            %% col1: hit / col2: corr_rej / col3: corr_rej-hit // row: subjects
            rewarding.(hpc_or_ctx{:}).(curr_names{r}){sbj_i,1}=mean(nonzeros(rewarding.hit.(hpc_or_ctx{:}){r,sbj_i}));
            rewarding.(hpc_or_ctx{:}).(curr_names{r}){sbj_i,2}=mean(nonzeros(rewarding.corr_rej.(hpc_or_ctx{:}){r,sbj_i}));
            
            rewarding.(hpc_or_ctx{:}).(curr_names{r}){sbj_i,3}=rewarding.(hpc_or_ctx{:}).(curr_names{r}){sbj_i,2}-rewarding.(hpc_or_ctx{:}).(curr_names{r}){sbj_i,1};
            
            rewarding.(hpc_or_ctx{:}).idx.corr_rej{sbj_i}=corr_rej_idx;
            rewarding.(hpc_or_ctx{:}).idx.hit{sbj_i}=hit_idx;
            rewarding.(hpc_or_ctx{:}).idx.false{sbj_i}=false_idx;
            rewarding.(hpc_or_ctx{:}).idx.miss{sbj_i}=miss_idx;
            
        end
    end
    % bilateral
    ori_name=strrep(bi_names,'Bi_','');
    for r=1:numel(bi_names)
        pat_L=cell2mat(rewarding.(hpc_or_ctx{:}).(curr_names{r}));
        pat_R=cell2mat(rewarding.(hpc_or_ctx{:}).(curr_names{r+(numel(curr_names)/2)}));
        
        pat_Bi=(pat_L+pat_R)/2;
        rewarding.(hpc_or_ctx{:}).(bi_names{r})=pat_Bi;
        
        
        % ROC로 나눠서 넣기
        for l=1:3
            lateral={'L_','R_','Bi_'};
            %                rewarding.contrast.EARLY.(hpc_or_ctx{:})=[curr_names, bi_names];
            %                rewarding.contrast.LATE.(hpc_or_ctx{:})=[curr_names, bi_names];
            
            pat={pat_L,pat_R,pat_Bi};
            rewarding.ROC_pattern.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(:,1);
            rewarding.ROC_pattern.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(:,2);
            rewarding.ROC_pattern.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(:,3);
            
            
            % performer group
            
            rewarding.ROC_pattern.EARLY.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,1);
            rewarding.ROC_pattern.EARLY.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,2);
            rewarding.ROC_pattern.EARLY.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,3);
            
            rewarding.ROC_pattern.LATE.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,1);
            rewarding.ROC_pattern.LATE.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,2);
            rewarding.ROC_pattern.LATE.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,3);
            
            rewarding.ROC_pattern.FAIL.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,1);
            rewarding.ROC_pattern.FAIL.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,2);
            rewarding.ROC_pattern.FAIL.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,3);


            rewarding.EARLY.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,1:2); %col 1: hit, col 2: corr_rej
            rewarding.LATE.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,1:2); %col 1: hit, col 2: corr_rej
            rewarding.FAIL.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,1:2); %col 1: hit, col 2: corr_rej

            rewarding.contrast.EARLY.(hpc_or_ctx{:})=[{strcat(lateral{l},ori_name{r})};num2cell(pat{l}(idx_early,3))];
            rewarding.contrast.LATE.(hpc_or_ctx{:})=[{strcat(lateral{l},ori_name{r})};num2cell(pat{l}(idx_late,3))];
            rewarding.contrast.FAIL.(hpc_or_ctx{:})=[{strcat(lateral{l},ori_name{r})};num2cell(pat{l}(idx_fail,3))];

        end
    end
    %%
    fn = fieldnames(rewarding.LATE.(hpc_or_ctx{:}));
    for fn_i = 1:length(fn)
        % EARLY contrast
        rewarding.contrast.EARLY.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.contrast.EARLY.(hpc_or_ctx{:})(2:length(idx_early)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.EARLY.subtract_corr_hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % LATE contrast
        rewarding.contrast.LATE.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.contrast.LATE.(hpc_or_ctx{:})(2:length(idx_late)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.LATE.subtract_corr_hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % FAIL contrast
        rewarding.contrast.FAIL.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.contrast.FAIL.(hpc_or_ctx{:})(2:length(idx_fail)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.FAIL.subtract_corr_hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % EARLY HIT
        rewarding.HIT.EARLY.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.HIT.EARLY.(hpc_or_ctx{:})(2:length(idx_early)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.EARLY.hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % LATE HIT
        rewarding.HIT.LATE.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.HIT.LATE.(hpc_or_ctx{:})(2:length(idx_late)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.LATE.hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % FAIL HIT
        rewarding.HIT.FAIL.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.HIT.FAIL.(hpc_or_ctx{:})(2:length(idx_fail)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.FAIL.hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % EARLY CORR_REJ
        rewarding.CORR_REJ.EARLY.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.CORR_REJ.EARLY.(hpc_or_ctx{:})(2:length(idx_early)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.EARLY.corr_rej.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % LATE CORR_REJ
        rewarding.CORR_REJ.LATE.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.CORR_REJ.LATE.(hpc_or_ctx{:})(2:length(idx_late)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.LATE.corr_rej.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
        % LATE CORR_REJ
        rewarding.CORR_REJ.FAIL.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
        rewarding.CORR_REJ.FAIL.(hpc_or_ctx{:})(2:length(idx_fail)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.ROC_pattern.FAIL.corr_rej.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
        
    end
    
end
save('half_2_corr+incorr_rewarding_patterns.mat',"rewarding")
%     save('half_2_rewarding_patterns.mat',"rewarding")
% end

%% excel
for perf={'EARLY','LATE','FAIL'}
    for reg={'hpc','ctx'}
        fn = fieldnames(rewarding.LATE.(reg{:}));
        for hcc={'HIT','CORR_REJ','contrast'}
            T = cell2table(rewarding.(hcc{:}).(perf{:}).(reg{:})(2:end,:),'VariableNames', fn');
            
            writetable(T, 'half_2_corr+incorr_rewarding_pattern.xlsx','Sheet',strcat(perf{:},'_',hcc{:},'_',reg{:}));
        end
    end
end

%% half 2-1

half_1=load('half_1_rewarding_patterns.mat');
half_2=load('half_2_rewarding_patterns.mat');

for perf={'EARLY','LATE','FAIL'}
    for reg={'hpc','ctx'}
        fn = fieldnames(rewarding.FAIL.(reg{:}));
        
        for hcc={'HIT','CORR_REJ','contrast'}
            hf1=cell2mat(half_1.rewarding.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
            hf2=cell2mat(half_2.rewarding.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
            
            T_cont = array2table(hf2-hf1,'VariableNames', fn');
            writetable(T_cont, 'HALF_rewarding_pattern.xlsx','Sheet',strcat('half2-1_', perf{:},'_',hcc{:},'_',reg{:}));
            
            T_half1 = array2table(hf1,'VariableNames', fn');
            writetable(T_half1, 'HALF_rewarding_pattern.xlsx','Sheet',strcat('HALF1_', perf{:},'_',hcc{:},'_',reg{:}));
            
            
            T_half2 = array2table(hf2,'VariableNames', fn');
            writetable(T, 'HALF_rewarding_pattern.xlsx','Sheet',strcat('HALF2_', perf{:},'_',hcc{:},'_',reg{:}));
            
        end
    end
end


%
%% same Ps, diff PS
% load('roi_name.mat')

for m_o = main_or_ODT
    disp(m_o)
    
    same_ctxt=struct;diff_ctxt=struct;

        load(string(fullfile(root_path,'spm_prep_glm_0922',m_o,'OC_pattern.mat')));

    
    % load('PS_basic_info.mat');
    % sbj_id_list=sbj_id_list(1:28)
    
    for hpc_or_ctx={'hpc','ctx'}
        if strcmp(hpc_or_ctx,'hpc')
            curr_names=roi_hpc_name(hpc_select);
        else
            curr_names=roi_ctx_name(ctx_select);
        end
        
        
        for sbj_i=1:numel(sbj_id_list)
            sbj_n = sbj_id_list(sbj_i);
            for r=1:numel(curr_names)
                
                curr_s_r=OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r};
                bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');
                
                if strcmp(m_o{:}, 'main')
                    
                    f_idx = find(cellfun(@(x) x == 1, curr_s_r(4,1:4)));
                    curr_f= cell2mat(OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,f_idx));
                    curr_c= setdiff(4:7,curr_f);
                    
                    f1_idx=find(cell2mat(curr_s_r(3,:))==curr_f(1));
                    f2_idx=find(cell2mat(curr_s_r(3,:))==curr_f(2));
                    c1_idx=find(cell2mat(curr_s_r(3,:))==curr_c(1));
                    c2_idx=find(cell2mat(curr_s_r(3,:))==curr_c(2));
                    
                else
                    f1_idx = find(cellfun(@(x) contains(x,'forest_1'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    f2_idx = find(cellfun(@(x) contains(x,'forest_2'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    c1_idx = find(cellfun(@(x) contains(x,'city_1'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    c2_idx = find(cellfun(@(x) contains(x,'city_2'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    
                end
                %% %%%%%%%%%%%%%%%%%%%%% 기존 within laps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% same
                %                 if issymmetric(bulk_corr)
                %                     tri_bulk=tril(bulk_corr,-1);
                %                 end
                
                %                 same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=[bulk_corr(f1_idx,f2_idx); bulk_corr(c1_idx,c2_idx)];
                beta_same=[bulk_corr(f1_idx,f2_idx); bulk_corr(c1_idx,c2_idx)];
                beta_same_re=reshape(beta_same,numel(beta_same),1);
                tri_bulk_corr=triu(bulk_corr,1);
                same_same= [tri_bulk_corr(f1_idx,f1_idx); tri_bulk_corr(c1_idx,c1_idx);tri_bulk_corr(f2_idx,f2_idx); tri_bulk_corr(c2_idx,c2_idx)];
                same_same=same_same(same_same~=0);
                if strcmp(m_o{:}, 'main')
                    %                     same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=beta_same;
                else
                    same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=[beta_same_re; same_same];
                end
                
                %                 if strcmp(m_o{:}, 'main')
                %                     m=[];
                %                     for lap=1:8
                %                         m(lap)=mean([same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}(lap,lap); same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}(lap+8,lap)],"all");
                %                     end
                %                     same_ctxt.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m';
                %                     for block=1:4
                %                         same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
                %                     end
                %                 end
                %% diff
                % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2
                
                %                 diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}=bulk_corr(f1_idx,c1_idx);
                %                 diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}=bulk_corr(f1_idx,c2_idx);
                %                 diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}=bulk_corr(c1_idx,f2_idx);
                %                 diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4}=bulk_corr(f2_idx,c2_idx);
                %
                
                %                 if strcmp(m_o{:}, 'main')
                %                     m=[];
                %                     for lap=1:8
                %                         % within lap
                %                         m(lap)=mean([diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}(lap,lap); diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}(lap,lap);diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}(lap,lap);diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4}(lap,lap)],"all");
                %                     end
                %                     diff_ctxt.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m';
                %                     for block=1:4
                %                         diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
                %                     end
                %                 end
                
                
                
                
                %% %%%%%%%%%%%%%%%%%%%%% across laps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(m_o{:}, 'main')
                    %% same
                    
                    bulk_same_f=bulk_corr(f1_idx,f2_idx);
                    bulk_same_c=bulk_corr(c1_idx,c2_idx);
                    same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=[bulk_same_f;bulk_same_c];
                    
                    m=[];
                    for lap=1:8
                        corr_f=bulk_same_f;
                        corr_c=bulk_same_c;
                        
                        within_f=corr_f(lap,lap);
                        within_c=corr_c(lap,lap);
                        
                        across_f=corr_f(lap,:);
                        across_c=corr_c(lap,:);
                        
                        across_f(across_f==within_f)=[];
                        across_c(across_c==within_c)=[];
                       
                        m_w(lap)=mean([within_f;within_c],"all");
                        m_a(lap)=mean([across_f;across_c],"all");
                    end
                    same_ctxt.within.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m_w';
                    same_ctxt.across.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m_a';
                    
                    for block=1:4
                        same_ctxt.within.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m_w(block*2-1:block*2));
                        same_ctxt.across.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m_a(block*2-1:block*2));
                    end
                    
                    %% diff
                    % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2
                    
                    diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}=bulk_corr(f1_idx,c1_idx);
                    diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}=bulk_corr(f1_idx,c2_idx);
                    diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}=bulk_corr(c1_idx,f2_idx);
                    diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4}=bulk_corr(f2_idx,c2_idx);
                    
                    
                    m=[];
                    for lap=1:8
                        
                        corr_f1=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1};
                        corr_f2=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2};
                        corr_c1=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3};
                        corr_c2=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4};
                        
                        within_f=[corr_f1(lap,lap);corr_f2(lap,lap)];
                        within_c=[corr_c1(lap,lap);corr_c2(lap,lap)];
                        
                        across_f=setdiff([corr_f1(lap,:);corr_f2(lap,:)],within_f);
                        across_c=setdiff([corr_c1(lap,:);corr_c2(lap,:)],within_c);
                        
                        m_w(lap)=mean([within_f;within_c],"all"); % within lap
                        m_a(lap)=mean([across_f;across_c],"all"); % across lap
                    end
                    diff_ctxt.within.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m_w';
                    diff_ctxt.across.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m_a';
                    
                    for block=1:4
                        diff_ctxt.within.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m_w(block*2-1:block*2));
                        diff_ctxt.across.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m_a(block*2-1:block*2));
                    end
                    
                    
                    
                end
            end
            
            
            if strcmp(m_o{:}, 'main')
                
                %% all_mean and half
                % within
                same_ctxt.within.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(same_ctxt.within.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                same_ctxt.within.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(same_ctxt.within.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
                diff_ctxt.within.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.within.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                diff_ctxt.within.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.within.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
                
                same_ctxt.within.(hpc_or_ctx{:}).all_mean{1,sbj_i}= mean(cell2mat(same_ctxt.within.(hpc_or_ctx{:}).half(1,sbj_i)));
                diff_ctxt.within.(hpc_or_ctx{:}).all_mean{1,sbj_i}= mean(cell2mat(diff_ctxt.within.(hpc_or_ctx{:}).half(1,sbj_i)));
                
                
                % across
                same_ctxt.across.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(same_ctxt.across.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                same_ctxt.across.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(same_ctxt.across.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
                diff_ctxt.across.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.across.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                diff_ctxt.across.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.across.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
                
                same_ctxt.across.(hpc_or_ctx{:}).all_mean{1,sbj_i}= mean(cell2mat(same_ctxt.across.(hpc_or_ctx{:}).half(1,sbj_i)));
                diff_ctxt.across.(hpc_or_ctx{:}).all_mean{1,sbj_i}= mean(cell2mat(diff_ctxt.across.(hpc_or_ctx{:}).half(1,sbj_i)));
                
            end
        end
        
        % save
     
            save(string(strcat(m_o{:},'_patterns_0922')),"same_ctxt","diff_ctxt","sbj_id_list_41")

        
    end
    
end
    %     %% across lap
    %
    %     for m_o = main_or_ODT
    %         % main_or_ODT = {'main'};
    %         disp(m_o{:});
    %
    %         if strcmp(m_o{:}, 'main')
    %             load(string(strcat(m_o,ver,'_patterns_0825')));
    %         else
    %             load(string(strcat(m_o,'_patterns_0825')));
    %         end
    %         for hpc_or_ctx={'hpc','ctx'}
    %             if strcmp(hpc_or_ctx,'hpc')
    %                 curr_names=roi_hpc_name(hpc_select);
    %             else
    %                 curr_names=roi_ctx_name(ctx_select);
    %             end
    %             for r=1:numel(curr_names)
    %                 for sbj_i=1:numel(sbj_id_list)
    %                     m2=[];
    %                     % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2
    %                     m2=mean([mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}); mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}); mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}); mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4})]);
    %                     diff_ctxt.new_mean.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m2;
    %
    %                     if strcmp(m_o{:}, 'main')
    %                         for block=1:4
    %                             diff_ctxt.new_mean.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m2(block*2-1:block*2));
    %                         end
    %                         diff_ctxt.new_mean.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
    %                         diff_ctxt.new_mean.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
    %                     end
    %                 end
    %             end
    %             for sbj_i=1:numel(sbj_id_list)
    %                 diff_ctxt.new_mean.(hpc_or_ctx{:}).all_mean{1,sbj_i}=mean(cell2mat(diff_ctxt.new_mean.(hpc_or_ctx{:}).half(1,sbj_i)));
    %             end
    %         end
    %         if strcmp(m_o{:}, 'main')
    %             save(string(strcat(m_o,ver,'_patterns_0825')),"same_ctxt","diff_ctxt","sbj_id_list_41")
    %         else
    %             save(string(strcat(m_o,'_patterns_0825')),"same_ctxt","diff_ctxt","sbj_id_list_41")
    %         end
    %     end
    
    %%
    
    now=struct;
    sbj_id_list=sbj_id_list_41;
    
    for m_o = main_or_ODT
        
        if strcmp(m_o{:}, 'main')
            load(string(strcat(m_o,'_patterns_0922')));
            p={'lap','half','block','all_mean'};
        else
            load(string(strcat(m_o,'_patterns_0922')));
            p=extractBefore(m_o{:},'-');
        end
        
        
        for diff_same={'diff','same'}
            for wi_ac={'within','across'}
                data={};
                for hpc_or_ctx={'hpc','ctx'}
                    if strcmp(hpc_or_ctx,'hpc')
                        curr_names=roi_hpc_name(hpc_select);
                    else
                        curr_names=roi_ctx_name(ctx_select);
                    end
                    for r=1:numel(curr_names)
                        for pp=p
                            curr_pp={};
                            
                            for sbj_i=1:numel(sbj_id_list)
                                sbj_n = sbj_id_list(sbj_i);
                                if strcmp(m_o{:}, 'main')
                                    temp_diff=diff_ctxt.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}){1,sbj_i}(:,r)';
                                    temp_same=same_ctxt.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}){1,sbj_i}(:,r)';
                                    %                                 temp_diff_across=diff_ctxt.new_mean.(hpc_or_ctx{:}).(pp{:}){1,sbj_i}(:,r)';
                                    % now.(hpc_or_ctx{:}).(pp{:}).same.(curr_names{r}){sbj_i,:}=temp_same;
                                    
                                    for ts=1:numel(temp_same)
                                        now.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}).same.(curr_names{r}){sbj_i,ts}=temp_same(ts);
                                        now.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}).diff.(curr_names{r}){sbj_i,ts}=temp_diff(ts);
                                        %                                     now.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}).diff_across.(curr_names{r}){sbj_i,ts}=temp_diff_across(ts);
                                        
                                        curr_pp{ts}=char(strcat(pp{:}, '_', string(ts)));
                                    end
                                    
                                    data = now.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}).(diff_same{:}).(curr_names{r});
                                    
                                    % now.diff.(curr_names{r}){sbj_i,:}=cell2mat(diff_ctxt.block  {1,sbj_i}(:,r)');
                                    % now.same.(curr_names{r}){sbj_i,:}=cell2mat(same_ctxt.block{1,sbj_i}(:,r)');
                                else
                                    now.(hpc_or_ctx{:}).(pp{:}).diff.(curr_names{r}){sbj_i,:}=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}),"all");
                                    now.(hpc_or_ctx{:}).(pp{:}).same.(curr_names{r}){sbj_i,:}=mean(same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i},"all");
                                    data(sbj_i, r) = now.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}).(diff_same{:}).(curr_names{r}){sbj_i};
                                end
                                
                            end
                            
                            if strcmp(m_o{:}, 'main')
                                T = cell2table(data, 'VariableNames', curr_pp);
                                writetable(T, string(strcat(m_o,'_',(wi_ac{:}),'_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet',curr_names{r});
                                
                            else
                                T = cell2table(data, 'VariableNames', curr_names);
                                writetable(T, string(strcat(m_o,'_final_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
                            end
                        end
                    end
                end
            end
        end
        save('now.mat','now')
    end 
   
        
        %% bilateral region
        
        for m_o = main_or_ODT
         for diff_same={'diff','same'}
            for wi_ac={'within','across'}
                data = [];
                for hpc_or_ctx={'hpc','ctx'}
                    if strcmp(hpc_or_ctx,'hpc')
                        curr_names=roi_hpc_name(hpc_select);
                    else
                        curr_names=roi_ctx_name(ctx_select);
                    end
                    bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
                    
                    if strcmp(m_o{:}, 'main')
                        p={'lap','half','block','all_mean'};
                        for pp=p
                            T_Bi=[];
                            
                            for r=1:numel(bi_names)
                                
                                T_L=readtable(string(strcat(m_o,'_',(wi_ac{:}),'_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet', curr_names{r});
                                T_R=readtable(string(strcat(m_o,'_',(wi_ac{:}),'_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet', curr_names{r+(numel(curr_names)/2)});
                                
                                
                                for col=1:numel(T_L(1,:))
                                    T_Bi(:,col)=mean([table2array(T_L(:,col)),table2array(T_R(:,col))],2);
                                end
                                
                                now.(wi_ac{:}).(hpc_or_ctx{:}).(pp{:}).(diff_same{:}).(bi_names{r})=T_Bi;
                                
                                T = array2table(T_Bi, 'VariableNames',T_L.Properties.VariableNames);
                                writetable(T, string(strcat(m_o,'_',(wi_ac{:}),'_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet',bi_names{r});
                                
                            end
                        end
                    else
                        T=readtable(string(strcat(m_o,'_final_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
                        
                        
                        for r=1:numel(bi_names)
                            data(:, r) = mean([table2array(T(:,r)),table2array(T(:,r+numel(bi_names)))],2);
                        end
                        
                        T_bi = array2table(data, 'VariableNames', bi_names);
                        T_combi = [T, T_bi];
                        
                        writetable(T_combi, string(strcat(m_o, '_final_pattern_0922.xlsx')), 'Sheet', strcat((diff_same{:}), '_', (hpc_or_ctx{:})));
                    end
                end
            end
            
        end
            save('now.mat','now')
        end
        
        %% if ODT, post-pre / if main, same-diff
        
        
        for hpc_or_ctx={'hpc','ctx'}
            if strcmp(hpc_or_ctx,'hpc')
                curr_names=roi_hpc_name(hpc_select);
            else
                curr_names=roi_ctx_name(ctx_select);
            end
            bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
            roi=[curr_names,bi_names];
            
            if ~strcmp(m_o{:},'main')
                
                for diff_same={'diff','same'}
                    
                    odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).subt=odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).post - odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).pre;
                    
                    T = array2table(odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).subt, 'VariableNames', roi);
                    
                    writetable(T, 'post-pre_final_pattern_0922.xlsx', 'Sheet', strcat((diff_same{:}), '_', (hpc_or_ctx{:})));
                end
                PSdiff=odt_ptn.(hpc_or_ctx{:}).same.subt - odt_ptn.(hpc_or_ctx{:}).diff.subt;
                T = array2table(PSdiff, 'VariableNames', roi);
                
                writetable(T, 'post-pre_final_pattern_0922.xlsx', 'Sheet', strcat('same-diff_', (hpc_or_ctx{:})));
                
                
            else
                p={'lap','half','block','all_mean'};
                for pp=p
                    
                    temp_tbl=string(strcat(m_o,'_',(wi_ac{:}),'_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx'));
                    sheets=sheetnames(temp_tbl);
                    for shts=1:numel(sheets)
                        same_ptn=readtable(string(strcat(m_o,'_',(wi_ac{:}),'_same_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet',sheets(shts));
                        diff_ptn=readtable(string(strcat(m_o,'_',(wi_ac{:}),'_diff_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet',sheets(shts));
                        
                        s_d_ptn=table2array(same_ptn)-table2array(diff_ptn);
                        
                        T = array2table(s_d_ptn, 'VariableNames', same_ptn.Properties.VariableNames);
                        writetable(T,string(strcat(m_o,'_',(wi_ac{:}),'_same-diff_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet',sheets(shts));
                    end
                end
            end
        end
        
        
        
        
        
        
        
           
        
        %% 3개 plot                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
%         load('main_acrossNwithin_patterns_0922');
        load('now.mat')
        curr_names={'Lt_Hp','Lt_DGCA3','Lt_CA1','Rt_Hp','Rt_DGCA3','Rt_CA1','Bi_Hp','Bi_DGCA3','Bi_CA1'};
        
        % data preparation
        %       p={'lap','half','block','all_mean'};wi_ac={'within','across'};
        
        
        for r=1:numel(curr_names)-3
            for sbj_i=1:numel(sbj_id_list)
                s_w=cell2mat(now.within.hpc.lap.same.(curr_names{r}){sbj_i,1});
                d_w=cell2mat(now.within.hpc.lap.diff.(curr_names{r}){sbj_i,1});
                s_a=cell2mat(now.across.hpc.lap.same.(curr_names{r}){sbj_i,1});
                d_a=cell2mat(now.across.hpc.lap.diff.(curr_names{r}){sbj_i,1});
                
                for l=1:8
                    same_within.(sprintf('lap_%d',l)){sbj_i,r}= double(s_w(l));
                    diff_within.(sprintf('lap_%d',l)){sbj_i,r}= double(d_w(l));
                    same_across.(sprintf('lap_%d',l)){sbj_i,r}= double(s_a(l));
                    diff_across.(sprintf('lap_%d',l)){sbj_i,r}= double(d_a(l));
                end
            end
        end

        for r=numel(curr_names)-2:numel(curr_names)
            for l=1:8
                for sbj_i=1:numel(sbj_id_list)
                    same_within.(sprintf('lap_%d',l)){sbj_i,r}= double(now.within.hpc.lap.same.(curr_names{r})(sbj_i,l));
                    diff_within.(sprintf('lap_%d',l)){sbj_i,r}= double(now.within.hpc.lap.diff.(curr_names{r})(sbj_i,l));
                    same_across.(sprintf('lap_%d',l)){sbj_i,r}= double(now.across.hpc.lap.same.(curr_names{r})(sbj_i,l));
                    diff_across.(sprintf('lap_%d',l)){sbj_i,r}= double(now.across.hpc.lap.diff.(curr_names{r})(sbj_i,l));
                end
            end
        end
        
        %%
        addpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code')
        
        curr_names= {'L.Hp','L.DGCA3','L.CA1','R.Hp','R.DGCA3','R.CA1','Bi.Hp','Bi.DGCA3','Bi.CA1'};
        % plotting!!
        for l=1:8
            within_diff=cell2mat(diff_within.(sprintf('lap_%d',l)));
            across_diff=cell2mat(diff_across.(sprintf('lap_%d',l)));
            within_same=cell2mat(same_within.(sprintf('lap_%d',l)));
            across_same=cell2mat(same_across.(sprintf('lap_%d',l)));
            
            M=[mean(within_diff); mean(across_diff); mean(within_same); mean(across_same)]';
            
            figure('Position',[928,291,1336,751]);
            ax = gca;
            h = bar(ax,M, 'grouped','LineWidth',1,'FaceColor','flat');
            hold(ax,'on');
            
            sem_w_d = std(within_diff) / sqrt(size(within_diff, 1));
            sem_a_d = std(across_diff) / sqrt(size(across_diff, 1));
            sem_w_s = std(within_same) / sqrt(size(within_same, 1));
            sem_a_s = std(across_same) / sqrt(size(across_same, 1));
            a=[sem_w_d; sem_a_d; sem_w_s; sem_a_s]';
            
            err(1)=errorbar(h(1).XEndPoints, mean(within_diff), sem_w_d, 'k', 'linestyle', 'none','LineWidth',1.5);
            err(2)=errorbar(h(2).XEndPoints, mean(across_diff), sem_a_d, 'k', 'linestyle', 'none','LineWidth',1.5);
            err(3)=errorbar(h(3).XEndPoints, mean(within_same), sem_w_s, 'k', 'linestyle', 'none','LineWidth',1.5);
            err(4)=errorbar(h(4).XEndPoints, mean(across_same), sem_a_s, 'k', 'linestyle', 'none','LineWidth',1.5);
            
            set(ax, 'XTick', 1:size(M,1), 'XTickLabel', curr_names,'LineWidth',1, 'FontSize',18,'FontWeight','bold');
            xtickangle(45); h(1).FaceColor = '#0072BD'; h(2).FaceColor = '#4DBEEE'; h(3).FaceColor = '#D95319';h(4).FaceColor = '#EDB120';
            legend(h, {'within diff', 'across diff', 'within same','across same'},'fontsize',9,'Location', 'bestoutside');
            box off; title(sprintf('lap %d diff vs same within vs across',l));
            
            %stats
            [~, p_w_d] = ttest(within_diff);
            [~, p_a_d] = ttest(across_diff);
            [~, p_w_s] = ttest(within_same);
            [~, p_a_s] = ttest(across_same);
            func_bar_significance(h, [p_w_d; p_a_d; p_w_s; p_a_s] , err)
            
            hold on;
            [~, p] = ttest(within_diff, across_diff);
            pval.diff.(sprintf('lap_%d',l))=p;
            idx=[]; idx=find(p<0.1); if ~isempty(idx); pval.signif.diff=curr_names(idx); end
            
            [~, p] = ttest(within_same, across_same);
            pval.same.(sprintf('lap_%d',l))=p;
            idx=[]; idx=find(p<0.1); if ~isempty(idx); pval.signif.same=curr_names(idx); end
            
            [~, p] = ttest(within_diff, across_diff);
            pval.within.(sprintf('lap_%d',l))=p;
            idx=[]; idx=find(p<0.1); if ~isempty(idx); pval.signif.within=curr_names(idx); end
            
            [~, p] = ttest(within_same, across_same);
            pval.across.(sprintf('lap_%d',l))=p;
            idx=[]; idx=find(p<0.1); if ~isempty(idx); pval.signif.across=curr_names(idx); end
            
            ylim([-0.15, 0.15]);
            hold off;
            
            
            saveas(gcf,['./across_within_comparison/' sprintf('lap_%d__diff_vs_same__within_vs_across.png',l)])
            
        end
