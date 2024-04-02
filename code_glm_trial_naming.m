%% set parameters
% load('regressors_GLM_0402')
sbj_id_list=1;
path_root = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\data_fmri_glm';


%% Input!
main_or_ODT = 'main'; % if you want to analize MAIN task, put 'main' or for ODT, put 'ODT'
threshold_type = 'FWE' ;% 'FWE' or 'FDR' or 'none'
date = '240330';
reg_mov_2or4 = 4; % if 4, moving regressor will contain [context(2)*order(first/remainder)], if 2, moving regressor will consider only context
phase = 'obj_show'; % input
if strcmp(main_or_ODT,'ODT'); phase ='ODT';end
if reg_mov_2or4==4; reg=reg_for_glm_ver1; elseif reg_mov_2or4==2; reg=reg_for_glm_ver2;end

addpath(genpath(fullfile(path_root,date)));

%% folder name
for i = 1: numel(sbj_id_list)
    n_sbj = sbj_id_list(i);
    c_sbj = sprintf('sub-%.2d',n_sbj);

    folders = dir(fullfile(path_root, date,(main_or_ODT),(phase),'1st_Level','single_trial',c_sbj,'betas','Sess001',strcat((phase),'*')));

    if strcmp(phase,'obj_show')
        suffix=[reg{1, i}.((main_or_ODT)).regress_suffix{1} reg{1, i}.((main_or_ODT)).regress_suffix{2}];
    elseif strcmp(phase,'moving')
        suffix=[reg{1, i}.((main_or_ODT)).regress_suffix{5} reg{1, i}.((main_or_ODT)).regress_suffix{6} reg{1, i}.((main_or_ODT)).regress_suffix{7} reg{1, i}.((main_or_ODT)).regress_suffix{8}];
    end

    %% copy and save folder
    addpath(genpath(folders(1).folder))
    cd(folders(1).folder)
    new_folder='../naming/'; mkdir(new_folder);

    for f=1:height(folders)
        oldname = folders(f).name;
        newname = strcat(folders(f).name,'_', suffix{f});
        status = copyfile(oldname, fullfile(new_folder, newname));
    end

    if status
        fprintf('copied folder successfully \n');
    else
        fprintf('fail copying folder\n');
    end

end



%% change beta file names

for sbj_i = 1: numel(sbj_id_list)
    n_sbj = sbj_id_list(sbj_i);
    c_sbj = sprintf('sub-%.2d',n_sbj);

    subfolders = dir(fullfile(path_root, date,(main_or_ODT),(phase),'1st_Level','single_trial',c_sbj,'betas','naming','*','*'));

    for ii = 1:length(subfolders)
        if subfolders(ii).isdir && ~strcmp(subfolders(ii).name,'..')

            load(fullfile(subfolders(ii).folder, subfolders(ii).name, 'SPM.mat'));
            betaFiles = dir(fullfile(subfolders(ii).folder, subfolders(ii).name, 'beta_*.nii'));

            for j = 1:length(betaFiles)
                betaNumber = sscanf(betaFiles(j).name, 'beta_%d.nii');
                if betaNumber >= 1 && betaNumber <= (length(betaFiles)-7)
                    newName = sprintf('beta_%04d_%s.nii', betaNumber, SPM.Sess.U(betaNumber).name{1});

                    if ~exist(fullfile(betaFiles(j).folder, newName), 'file')

                        movefile(fullfile(betaFiles(j).folder, betaFiles(j).name), fullfile(betaFiles(j).folder, newName));
                    end
                end
            end
        end

    end
end














%
%
%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % v1_onsets 생성
% v1_onsets = {sort([onsets{1, 1} onsets{1, 4}]); sort([onsets{1, 2} onsets{1, 3} onsets{1, 5}])};
%
% % 파일 이름 설정
% output_filename = 'output.txt';
%
% % 파일 열기
% fileID = fopen(output_filename, 'w');
%
% % 각각의 v1_onsets에서 hit 또는 corr_rej인지 확인하고 인덱스 추가하여 파일에 저장
% for i = 1:length(v1_onsets)
%     for j = 1:length(v1_onsets{i})
%         onset = v1_onsets{i}(j);
%         if ismember(onset, onsets{1, 1}) || ismember(onset, onsets{1, 4})
%             fprintf(fileID, 'beta_%04d_hit\n', onset);
%         else
%             fprintf(fileID, 'beta_%04d_corr_rej\n', onset);
%         end
%     end
% end
%
% % 파일 닫기
% fclose(fileID);
%
% disp(['인덱스가 추가된 파일이 ' output_filename '에 저장되었습니다.']);
%
%
%
%
%
