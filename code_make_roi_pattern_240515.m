clear;clc;

cd 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code'
addpath('C:\Users\User\Documents\MATLAB');
load('sbj_id_list');


TYPE='FS60'; % 'HBT' or 'FS60' or 'extra'
make_pattern=1; % if 1, pattern will be made
do_analysis=1; % if u want to do pattern analysis, type 1
mask_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_bids\derivatives\HPC_seg\' TYPE];
glm_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_glm\240529\ODT\1st_Level';
out_path=fullfile(glm_path,'../pattern_240529');mkdir(out_path);
addpath(genpath(glm_path))

%% make pattern - ver1
% roi=struct;
% for sbj_i=1:numel(sbj_id_list)
%     sbj_n = sbj_id_list(sbj_i);
%     local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'coreg_*.nii'));
%
%     beta_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_glm\240426_anat\ODT\1st_Level\' sprintf('sub-%.2d',sbj_n)];
%     betas=dir(fullfile(beta_path,'coreg_beta_*.nii'));
%
%     for b=1:height(betas)
%         curr_beta = fullfile(betas(1).folder,betas(b).name);
%         for m=1:height(local_masks)
%             curr_mask=niftiread(fullfile(local_masks(1).folder,local_masks(m).name));
%             roi_coord = mask2coord(curr_mask);
%
%             roi_beta = niftiread(curr_beta);
%             roi_beta_2d = reshape(roi_beta,416,354);% roi_beta reshape할 때 beta 마다 달라짐!!!
%             %ODT: 52*59*48 = 416*354 = 147264
%             %obj_show: 1001*300 = 300300
%
%             roi.(sprintf('s%.2d',sbj_n)).hpc.(TYPE).roi_name_list{m}=extractAfter(local_masks(m).name,'coreg_');
%             roi.(sprintf('s%.2d',sbj_n)).hpc.(TYPE).roi_pattern{m}=roi_beta_2d;
%         end
%     end
% end


% %% make pattern - ver2
% if make_pattern==1
%     ODT_pattern=struct;
%     for sbj_i=1:numel(sbj_id_list)
%         sbj_n = sbj_id_list(sbj_i);
%
%         local_masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'coreg_*.nii'));
%
%         beta_path=['Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_glm\240426_anat\ODT\1st_Level\' sprintf('sub-%.2d',sbj_n)];
%         betas=dir(fullfile(beta_path,'coreg_beta_*.nii'));
%
%
%         for m=1:height(local_masks)
%             curr_mask=niftiread(fullfile(local_masks(1).folder,local_masks(m).name));
%             roi_coord = mask2coord(curr_mask);
%             for b=1:height(betas)-7
%                 curr_beta = fullfile(betas(1).folder,betas(b).name);
%                 % beta_name
%                 load(fullfile(betas(1).folder,'SPM.mat'))
%                 curr_beta_name = extractBetween(SPM.Vbeta(b).descrip,"Sn(1) ","*bf(1)");
%
%                 roi_beta = niftiread(curr_beta);
%                 roi_beta_2d = reshape(roi_beta,416,354);% roi_beta reshape할 때 beta 마다 달라짐!!!
%                 %ODT: 52*59*48 = 416*354 = 147264
%                 %obj_show: 1001*300 = 300300
%
%                 roi_name = extractBetween(local_masks(m).name,'coreg_','.nii');
%                 ODT_pattern.(TYPE).roi_name_list{m}=roi_name{1};
%                 ODT_pattern.beta_name_list{b}=curr_beta_name{1};
%                 ODT_pattern.(TYPE).(sprintf('s%.2d',sbj_n)).(roi_name{1}){1,b}=curr_beta_name;
%                 ODT_pattern.(TYPE).(sprintf('s%.2d',sbj_n)).(roi_name{1}){2,b}=roi_beta_2d;
%             end
%         end
%     end
%     save("ODT_pattern");
% end

%% make pattern - ver3
%
% if make_pattern==1
%
% beta_names=dir(fullfile(glm_path,'roi_mask','extra'));
% beta_names=arrayfun(@(x) beta_names(x).name, 3:height(beta_names), 'UniformOutput', false)';
%
%
% ODT_pattern=struct;roi_names=struct;
% for sbj_i=1:numel(sbj_id_list)
%     sbj_n = sbj_id_list(sbj_i);
%     roi_name=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'*.nii'));
%     roi_name=arrayfun(@(x) extractBefore(roi_name(x).name,'.nii'), 1:height(roi_name), 'UniformOutput', false)';
%
%
%     for b=1:numel(beta_names)
%         beta_masks_path=fullfile(glm_path,'roi_mask', TYPE, beta_names{b},sprintf('sub-%.2d',sbj_n));
%
%         betas=dir(fullfile(glm_path,sprintf('sub-%.2d',sbj_n),'beta_*.nii'));
%         curr_beta = niftiread(fullfile(betas(1).folder,betas(b).name));
%
%         for r=1:numel(roi_name)
%             curr_roi_mask=niftiread(fullfile(beta_masks_path,strcat(roi_name{r},'.nii')));
%             pattern_3d=curr_beta;
%             pattern_3d(curr_roi_mask==0)=NaN;
%
%             % 2d로 변환
%             sz=size(pattern_3d);
%             pattern_2d=reshape(pattern_3d,[sz(1)*sz(2),sz(3)]);
%
%             %1d로 변환
%             pattern_1d=reshape(pattern_3d,[],1);
%
%              ODT_pattern.(TYPE).(beta_names{b}).(roi_name{r}){sbj_i}=pattern_2d;
%              ODT_pattern.(TYPE).pattern_3d.(beta_names{b}).(roi_name{r}){sbj_i}=pattern_3d;
%             roi_names.(TYPE){r,1}=roi_name{r};
%         end
%     end
% end
%     save(fullfile(glm_path,'../raw_ODT_pattern_240515'),"ODT_pattern","roi_names","beta_names");
% end

%%
% %% make pattern - ver4
%
% if make_pattern==1
%
%     beta_names=dir(fullfile(glm_path,'roi_mask_new',TYPE));
%     beta_names=arrayfun(@(x) beta_names(x).name, 3:height(beta_names), 'UniformOutput', false)';
%
%
%     ODT_pattern=struct;roi_names=struct;
%     for sbj_i=1:numel(sbj_id_list)
%         sbj_n = sbj_id_list(sbj_i);
%         roi_name=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'*.nii'));
%         roi_name=arrayfun(@(x) extractBefore(roi_name(x).name,'.nii'), 1:height(roi_name), 'UniformOutput', false)';
%
%
%         for b=1:numel(beta_names)
%             beta_masks_path=fullfile(glm_path,'roi_mask_new', TYPE, beta_names{b},sprintf('sub-%.2d',sbj_n));
%
%             betas=dir(fullfile(glm_path,sprintf('sub-%.2d',sbj_n),'beta_*.nii'));
%             curr_beta = niftiread(fullfile(betas(1).folder,betas(b).name));
%
%             for r=1:numel(roi_name)
%                 curr_roi_mask=niftiread(fullfile(beta_masks_path,strcat(beta_names{b},'_', roi_name{r},'.nii')));
%                 pattern_roi=curr_beta(curr_roi_mask==1);
%                 roiname=strrep(roi_name{r},'GC-','');
%                 ODT_pattern.(TYPE).(beta_names{b}).(roiname){sbj_i}=pattern_roi;
%                 roi_names.(TYPE){r,1}=roi_name{r};
%             end
%         end
%     end
%     save(fullfile(out_path,'raw_ODT_pattern'),"ODT_pattern","roi_names","beta_names");
% end
%

%%
%% make pattern - ver5

if make_pattern==1
load(fullfile(glm_path,'../names.mat'))
    beta_names=dir(fullfile(glm_path,'roi_mask_new',TYPE));
    beta_names=arrayfun(@(x) beta_names(x).name, 3:height(beta_names), 'UniformOutput', false)';


    ODT_pattern=struct;roi_names=struct;
    for sbj_i=2:numel(sbj_id_list)
        sbj_n = sbj_id_list(sbj_i);
        masks=dir(fullfile(mask_path,sprintf('s-%.2d',sbj_n),'*.nii'));
        roi_name=arrayfun(@(x) extractBefore(masks(x).name,'.nii'), 1:height(masks), 'UniformOutput', false)';

        for b=1:numel(beta_names)
            beta_path=fullfile(glm_path,'roi_mask_new', TYPE, beta_names{b},sprintf('sub-%.2d',sbj_n));

            for r=1:numel(roi_name)
                find_beta=dir(fullfile(beta_path,strcat('*_', roi_name{r},'.nii')));
                curr_beta=niftiread(fullfile(find_beta(1).folder,find_beta(1).name));
                curr_roi_mask=niftiread(fullfile(mask_path,sprintf('s-%.2d',sbj_n), strcat(roi_name{r},'.nii')));

                pattern_roi=curr_beta(curr_roi_mask);

                roiname=strrep(roi_name{r},'GC-','');
                ODT_pattern.(TYPE).(beta_names{b}).(roiname){sbj_i}=pattern_roi;
                roi_names.(TYPE){r,1}=roi_name{r};
            end
        end
    end
    save(fullfile(out_path,'raw_ODT_pattern'),"ODT_pattern","roi_names","beta_names");
end



%% correlation coefficient - ver2


% {'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2','post_forest_1','post_forest_2','post_city_1','post_city_2'}

%0. obje 4개 pre vs post ['f1','f2','c1','c2']
%1. city - pre vs post [-'pre_c1c2' +'post_c1c2']
%2. forest - pre vs post [-'pre_f1f2' +'postf1f2']
%3. forest vs city [('f1'+'f2')-('c1'+'c2')]
%4. pre - forest vs city ['pre_f1f2'-'pre_c1c2']
%5. post - forest vs city ['postf1f2'-'post_c1c2']


%'pre_forest_1''post_forest_1' + %'pre_forest_2' 'post_forest_2'  [forest] vs [city]

%('pre_forest_1''pre_forest_2')-('post_forest_1'post_forest_2' ) forest - [pre vs post]
%-(('pre_forest_1''pre_forest_2')('pre_forest_1''pre_forest_2')) +(('post_forest_1'post_forest_2' )('post_forest_1'post_forest_2' )) [pre vs post]

if do_analysis==1

    load(fullfile(out_path,'raw_ODT_pattern'))
    corr=struct;
    corr_name_list={'f1','f2','c1','c2','pre_f1f2','pre_c1c2','post_f1f2','post_c1c2','pre_f1c1','pre_f1c2','pre_f2c1','pre_f2c2','post_f1c1','post_f1c2','post_f2c1','post_f2c2'};
    bn={'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2','post_forest_1','post_forest_2','post_city_1','post_city_2'};

    for sbj_i=1:numel(sbj_id_list)
        sbj_n = sbj_id_list(sbj_i);
        for r=1:numel(roi_names.(TYPE))
            rn= roi_names.(TYPE){r};
            rn= strrep(rn,'GC-','');
            ptrn=ODT_pattern.(TYPE);
            ODT_pattern.(TYPE).corr.corr_name_list=corr_name_list;
            %% corrcoef
            R=corrcoef(ptrn.(bn{1}).(rn){sbj_i},ptrn.(bn{5}).(rn){sbj_i}, 'rows','pairwise');
            corr.f1=R(2);
            R=corrcoef(ptrn.(bn{2}).(rn){sbj_i},ptrn.(bn{6}).(rn){sbj_i}, 'rows','pairwise');
            corr.f2=R(2);
            R=corrcoef(ptrn.(bn{3}).(rn){sbj_i},ptrn.(bn{7}).(rn){sbj_i}, 'rows','pairwise');
            corr.c1=R(2);
            R=corrcoef(ptrn.(bn{4}).(rn){sbj_i},ptrn.(bn{8}).(rn){sbj_i}, 'rows','pairwise');
            corr.c2=R(2);
            R=corrcoef(ptrn.(bn{1}).(rn){sbj_i},ptrn.(bn{2}).(rn){sbj_i}, 'rows','pairwise');
            corr.pre_f1f2=R(2);
            R=corrcoef(ptrn.(bn{3}).(rn){sbj_i},ptrn.(bn{4}).(rn){sbj_i}, 'rows','pairwise');
            corr.pre_c1c2=R(2);
            R=corrcoef(ptrn.(bn{5}).(rn){sbj_i},ptrn.(bn{6}).(rn){sbj_i}, 'rows','pairwise');
            corr.post_f1f2=R(2);
            R=corrcoef(ptrn.(bn{7}).(rn){sbj_i},ptrn.(bn{8}).(rn){sbj_i}, 'rows','pairwise');
            corr.post_c1c2=R(2);

            %%
            R=corrcoef(ptrn.(bn{1}).(rn){sbj_i},ptrn.(bn{3}).(rn){sbj_i}, 'rows','pairwise');
            corr.pre_f1c1=R(2);
            R=corrcoef(ptrn.(bn{1}).(rn){sbj_i},ptrn.(bn{4}).(rn){sbj_i}, 'rows','pairwise');
            corr.pre_f1c2=R(2);
            R=corrcoef(ptrn.(bn{2}).(rn){sbj_i},ptrn.(bn{3}).(rn){sbj_i}, 'rows','pairwise');
            corr.pre_f2c1=R(2);
            R=corrcoef(ptrn.(bn{2}).(rn){sbj_i},ptrn.(bn{4}).(rn){sbj_i}, 'rows','pairwise');
            corr.pre_f2c2=R(2);
            R=corrcoef(ptrn.(bn{5}).(rn){sbj_i},ptrn.(bn{7}).(rn){sbj_i}, 'rows','pairwise');
            corr.post_f1c1=R(2);
            R=corrcoef(ptrn.(bn{5}).(rn){sbj_i},ptrn.(bn{8}).(rn){sbj_i}, 'rows','pairwise');
            corr.post_f1c2=R(2);
            R=corrcoef(ptrn.(bn{6}).(rn){sbj_i},ptrn.(bn{7}).(rn){sbj_i}, 'rows','pairwise');
            corr.post_f2c1=R(2);
            R=corrcoef(ptrn.(bn{6}).(rn){sbj_i},ptrn.(bn{8}).(rn){sbj_i}, 'rows','pairwise');
            corr.post_f2c2=R(2);


            %% target
            R=corrcoef(ptrn.pre_target.(rn){sbj_i},ptrn.post_target.(rn){sbj_i}, 'rows','pairwise');
            corr.target=R(2);



            for c = 1:numel(corr_name_list)
                ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(rn){1,c}=corr_name_list{c};
                ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(rn){2,c}=getfield(corr, corr_name_list{c});
            end

        end
    end
    save(fullfile(out_path,'ODT_pattern'),"ODT_pattern","roi_names","beta_names");

    %%
    % %% contrast - ver1
    % ODT_pattern.contrast.contrast_list={("0. obj 4개 pre vs post ['f1','f2','c1','c2']");("1. city - pre vs post [-'pre_c1c2' +'post_c1c2']");("2. forest - pre vs post [-'pre_f1f2' +'postf1f2']");...
    %     ("3. forest vs city [('f1'+'f2')-('c1'+'c2')]");("4. pre - forest vs city ['pre_f1f2'-'pre_c1c2']");("5. post - forest vs city ['postf1f2'-'post_c1c2']")};
    %
    % for sbj_i=1:numel(sbj_id_list)
    %     sbj_n = sbj_id_list(sbj_i);
    %     for r=1:numel(roi_names.(TYPE))
    %         rn= roi_names.(TYPE){r};
    %         corr_value=ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(rn)(2,:);
    %         ODT_pattern.contrast.(rn).pre_same(sbj_i) = mean([corr_value{5}, corr_value{6}]);
    %         ODT_pattern.contrast.(rn).post_same(sbj_i) = mean([corr_value{7}, corr_value{8}]);
    %         ODT_pattern.contrast.(rn).pre_diff(sbj_i) = mean([corr_value{9}, corr_value{10},corr_value{11}, corr_value{12}]);
    %         ODT_pattern.contrast.(rn).post_diff(sbj_i) = mean([corr_value{13}, corr_value{14},corr_value{15}, corr_value{16}]);
    %
    %         writematrix(ODT_pattern.contrast.(rn).pre_same', fullfile(glm_path,'../pattern','pre_same.xlsx'), 'Sheet', rn);
    %         writematrix(ODT_pattern.contrast.(rn).post_same', fullfile(glm_path,'../pattern','post_same.xlsx'), 'Sheet', rn);
    %         writematrix(ODT_pattern.contrast.(rn).pre_diff', fullfile(glm_path,'../pattern','pre_diff.xlsx'), 'Sheet', rn);
    %         writematrix(ODT_pattern.contrast.(rn).post_diff', fullfile(glm_path,'../pattern','post_diff.xlsx'), 'Sheet', rn);
    %
    %     end
    % end
    %
    % save(fullfile(glm_path,'../ODT_pattern_240515'),"ODT_pattern","roi_names","beta_names");


    %% contrast - ver2
    ODT_pattern.contrast.contrast_list={("0. obj 4개 pre vs post ['f1','f2','c1','c2']");("1. city - pre vs post [-'pre_c1c2' +'post_c1c2']");("2. forest - pre vs post [-'pre_f1f2' +'postf1f2']");...
        ("3. forest vs city [('f1'+'f2')-('c1'+'c2')]");("4. pre - forest vs city ['pre_f1f2'-'pre_c1c2']");("5. post - forest vs city ['postf1f2'-'post_c1c2']")};

    for r=1:numel(roi_names.(TYPE))
        rn= roi_names.(TYPE){r};
        rn= strrep(rn,'GC-','');

        M = zeros(numel(sbj_id_list), 4);

        for sbj_i=1:numel(sbj_id_list)
            sbj_n = sbj_id_list(sbj_i);
            corr_value=ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(rn)(2,:);

            M(sbj_i, 1) = mean([corr_value{5}, corr_value{6}]);  % pre_same
            M(sbj_i, 2) = mean([corr_value{7}, corr_value{8}]);  % post_same
            M(sbj_i, 3) = mean([corr_value{9}, corr_value{10},corr_value{11}, corr_value{12}]);  % pre_diff
            M(sbj_i, 4) = mean([corr_value{13}, corr_value{14},corr_value{15}, corr_value{16}]);  % post_diff
        end
        M_cell = [{'pre_same', 'post_same', 'pre_diff', 'post_diff'}; num2cell(M)];

        writecell(M_cell, fullfile(glm_path,'../pattern','contrast_pattern.xlsx'), 'Sheet', rn);
    end

    save(fullfile(out_path,'ODT_pattern'),"ODT_pattern","roi_names","beta_names");






    %% make excel file

    T=[];
    for c = 1:numel(corr_name_list)
        T = array2table(zeros(numel(sbj_id_list), numel(roi_names.(TYPE))), 'VariableNames', roi_names.(TYPE));

        for r=1:numel(roi_names.(TYPE))
            rn= roi_names.(TYPE){r};
            for sbj_i=1:numel(sbj_id_list)
                sbj_n = sbj_id_list(sbj_i);

                T.(rn)(sbj_i) = ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(rn){2,c};
            end
        end

        C = [roi_names.(TYPE)' ; table2cell(T)];
        writecell(C, fullfile(out_path, [TYPE '_patterns.xlsx']), 'Sheet', corr_name_list{c});
    end


end




















% %% correlation coefficient - ver1
%
%
% % {'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2','post_forest_1','post_forest_2','post_city_1','post_city_2'}
%
% %0. obje 4개 pre vs post ['f1','f2','c1','c2']
% %1. city - pre vs post [-'pre_c1c2' +'post_c1c2']
% %2. forest - pre vs post [-'pre_f1f2' +'postf1f2']
% %3. forest vs city [('f1'+'f2')-('c1'+'c2')]
% %4. pre - forest vs city ['pre_f1f2'-'pre_c1c2']
% %5. post - forest vs city ['postf1f2'-'post_c1c2']
%
%
% %'pre_forest_1''post_forest_1' + %'pre_forest_2' 'post_forest_2'  [forest] vs [city]
%
% %('pre_forest_1''pre_forest_2')-('post_forest_1'post_forest_2' ) forest - [pre vs post]
% %-(('pre_forest_1''pre_forest_2')('pre_forest_1''pre_forest_2')) +(('post_forest_1'post_forest_2' )('post_forest_1'post_forest_2' )) [pre vs post]
%
% if do_analysis==1
%
%     load('ODT_pattern')
%     corr=struct;
%     corr_name_list={'f1','f2','c1','c2','pre_f1f2','pre_c1c2','post_f1f2','post_c1c2','pre_f1c1','pre_f1c2','pre_f2c1','pre_f2c2','post_f1c1','post_f1c2','post_f2c1','post_f2c2'};
%
%     for sbj_i=1:numel(sbj_id_list)
%         sbj_n = sbj_id_list(sbj_i);
%         for m=1:numel(ODT_pattern.(TYPE).roi_name_list)
%             roi_name=ODT_pattern.(TYPE).roi_name_list{m};
%             betas = ODT_pattern.(TYPE).(sprintf('s%.2d',sbj_n)).(roi_name)(2,:);
%
%             ODT_pattern.corr.corr_name_list=corr_name_list;
%             %% corrcoef
%             R=corrcoef(betas{1},betas{5}, 'rows','pairwise');
%             corr.f1=R(2);
%             R=corrcoef(betas{2},betas{6}, 'rows','pairwise');
%             corr.f2=R(2);
%             R=corrcoef(betas{3},betas{7}, 'rows','pairwise');
%             corr.c1=R(2);
%             R=corrcoef(betas{4},betas{8}, 'rows','pairwise');
%             corr.c2=R(2);
%             R=corrcoef(betas{1},betas{2}, 'rows','pairwise');
%             corr.pre_f1f2=R(2);
%             R=corrcoef(betas{3},betas{4}, 'rows','pairwise');
%             corr.pre_c1c2=R(2);
%             R=corrcoef(betas{5},betas{6}, 'rows','pairwise');
%             corr.post_f1f2=R(2);
%             R=corrcoef(betas{7},betas{8}, 'rows','pairwise');
%             corr.post_c1c2=R(2);
%
%             %%
%             R=corrcoef(betas{1},betas{3}, 'rows','pairwise');
%             corr.pre_f1c1=R(2);
%             R=corrcoef(betas{1},betas{4}, 'rows','pairwise');
%             corr.pre_f1c2=R(2);
%             R=corrcoef(betas{2},betas{3}, 'rows','pairwise');
%             corr.pre_f2c1=R(2);
%             R=corrcoef(betas{2},betas{4}, 'rows','pairwise');
%             corr.pre_f2c2=R(2);
%             R=corrcoef(betas{5},betas{7}, 'rows','pairwise');
%             corr.post_f1c1=R(2);
%             R=corrcoef(betas{5},betas{8}, 'rows','pairwise');
%             corr.post_f1c2=R(2);
%             R=corrcoef(betas{6},betas{7}, 'rows','pairwise');
%             corr.post_f2c1=R(2);
%             R=corrcoef(betas{6},betas{8}, 'rows','pairwise');
%             corr.post_f2c2=R(2);
%
%             for c = 1:numel(corr_name_list)
%                 ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(roi_name){1,c}=corr_name_list{c};
%                 ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(roi_name){2,c}=getfield(corr, corr_name_list{c});
%             end
%
%         end
%     end
%     % save("ODT_pattern");
%     pattern_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\ocat_fmri_bids\derivatives\HPC_subregion\pattern';
%     save(fullfile(pattern_path,"ODT_pattern"));
%
%
%     %% contrast
%     ODT_pattern.corr.contrast_list={("0. obj 4개 pre vs post ['f1','f2','c1','c2']");("1. city - pre vs post [-'pre_c1c2' +'post_c1c2']");("2. forest - pre vs post [-'pre_f1f2' +'postf1f2']");...
%         ("3. forest vs city [('f1'+'f2')-('c1'+'c2')]");("4. pre - forest vs city ['pre_f1f2'-'pre_c1c2']");("5. post - forest vs city ['postf1f2'-'post_c1c2']")};
%
%     for sbj_i=1:numel(sbj_id_list)
%         sbj_n = sbj_id_list(sbj_i);
%         for m=1:numel(ODT_pattern.(TYPE).roi_name_list)
%             roi_name=ODT_pattern.(TYPE).roi_name_list{m};
%             corr_name=ODT_pattern.corr.(sprintf('s%.2d',sbj_n)).(roi_name)(2,:);
%             ODT_pattern.corr.(roi_name).pre_same(sbj_i) = mean([corr_name{5}, corr_name{6}]);
%             ODT_pattern.corr.(roi_name).post_same(sbj_i) = mean([corr_name{7}, corr_name{8}]);
%             ODT_pattern.corr.(roi_name).pre_diff(sbj_i) = mean([corr_name{9}, corr_name{10},corr_name{11}, corr_name{12}]);
%             ODT_pattern.corr.(roi_name).post_diff(sbj_i) = mean([corr_name{13}, corr_name{14},corr_name{15}, corr_name{16}]);
%
%             writematrix(ODT_pattern.corr.(roi_name).pre_same, fullfile(pattern_path,'pre_same.xlsx'), 'Sheet', roi_name);
%             writematrix(ODT_pattern.corr.(roi_name).post_same, fullfile(pattern_path,'post_same.xlsx'), 'Sheet', roi_name);
%             writematrix(ODT_pattern.corr.(roi_name).pre_diff, fullfile(pattern_path,'pre_diff.xlsx'), 'Sheet', roi_name);
%             writematrix(ODT_pattern.corr.(roi_name).post_diff, fullfile(pattern_path,'post_diff.xlsx'), 'Sheet', roi_name);
%
%
%
%             %
%             %         T_pre_same = array2table(ODT_pattern.corr.(roi_name).pre_same);
%             %         T_post_same = array2table(ODT_pattern.corr.(roi_name).post_same);
%             %         T_pre_diff = array2table(ODT_pattern.corr.(roi_name).pre_diff);
%             %         T_post_diff = array2table(ODT_pattern.corr.(roi_name).post_diff);
%             %
%             %         writetable(T_pre_same, 'pre_same.xlsx', 'Sheet', roi_name);
%             %         writetable(T_post_same, 'post_same.xlsx', 'Sheet', roi_name);
%             %         writetable(T_pre_diff, 'pre_diff.xlsx', 'Sheet', roi_name);
%             %         writetable(T_post_diff, 'post_diff.xlsx', 'Sheet', roi_name);
%
%         end
%     end
%
%     save(fullfile(pattern_path,"ODT_pattern"));
%
% end
%
%
%
%
