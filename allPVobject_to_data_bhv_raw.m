% 폴더 경로 설정
folder_path = 'C:\Users\Leelab\Downloads\drive-download-20231127T052723Z-001\';
subfolders = dir(fullfile(folder_path, 'Sub*'));
new_folder_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_pvobj\';

addpath(folder_path,new_folder_path);

% 각 Sub 폴더에 대해
for i = 1:length(subfolders)
    % pre-PVtaskLog.xlsx과 post-PVtaskLog.xlsx 파일 로드
    pre_file = readtable(fullfile(folder_path, subfolders(i).name, [subfolders(i).name '_pre-PVtaskLog.xlsx']));
    post_file = readtable(fullfile(folder_path, subfolders(i).name, [subfolders(i).name '_post-PVtaskLog.xlsx']));
    
    % 두 파일을 세로로 연결
    combined_file = [pre_file; post_file];
    
    % 새로운 .xlsx 파일로 저장
    writetable(combined_file, fullfile(new_folder_path, ['sub-' num2str(i, '%02.f') '_combined.csv']));
end


% data_bhv_raw 폴더 경로 설정
data_bhv_raw_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_raw\';
pv_data_bhv_raw_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\data_bhv_raw_pv\';
addpath(data_bhv_raw_path,pv_data_bhv_raw_path);

% 각 Sub 폴더에 대해
for i = 1:length(subfolders)
    % combined.csv 파일 로드
    combined_file = readtable(fullfile(new_folder_path, ['sub-' num2str(i, '%02.f') '_combined.csv']));
    
    % data_bhv_raw의 i번째 파일 로드
    fileToRead1= dir(fullfile(data_bhv_raw_path, '*Time_Behavior.csv'));
    data_bhv_raw_file = readtable(fileToRead1(i).name);
    

    % combined_file.Var2의 모든 값에 value_to_add 더하기
    first_MR_row = find(strcmp(data_bhv_raw_file.Var1, 'MR'), 1, 'first');
    value_to_add = data_bhv_raw_file.Var2(first_MR_row);
    combined_file.Var2 = combined_file.Var2 + value_to_add;


    % Var1 값이 같은 행을 찾습니다.
    idx = find(contains(combined_file.Var1,"Button"));
combined_file(idx,:)=[];
    idx = find(contains(data_bhv_raw_file.Var1,combined_file.Var1));
    for j = 1:length(idx)
    % Var4 값을 변경합니다.
    data_bhv_raw_file.Var4(idx(j))= combined_file.Var4(j);
   
    end


    % 수정된 data_bhv_raw 파일 저장
    writetable(data_bhv_raw_file, fullfile(pv_data_bhv_raw_path, fileToRead1(i).name));
end




