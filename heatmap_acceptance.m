root_path = 'C:\Users\leelab\Desktop\fMRI\lecture\ds000005';
task_name = 'mixedgamblestask';
tmp_folders = dir(fullfile('C:\Users\leelab\Desktop\fMRI\lecture\ds000005','sub-*')); 
dfolders = tmp_folders([tmp_folders(:).isdir]);
subj_list = {dfolders.name};

all_data = [];

for i = 1:length(subj_list)
    subjID = subj_list{i};

    func_path = fullfile(root_path, subjID, 'ses-01', 'func'); %%% CHECK if this path is correct

    % path for confounding factors
    event_file = dir(fullfile(func_path, ['*task-', task_name, '*events.tsv'])); 

    for j = 1:length(event_file)
        % 파일 이름 가져오기
        file_name = fullfile(event_file(j).folder, event_file(j).name);

        % 파일 읽기
        T = readtable(file_name, 'FileType', 'text', 'Delimiter', '\t');

        % 6, 7, 9열만 선택
        data = [T(:,6), T(:,7), T(:,9)];

        % all_data에 추가
        all_data = [all_data; data];
    end
end

% x축과 y축 범위 설정
x_range = 10:2:40;
y_range = 20:-1:5;

% 히트맵을 위한 빈 셀 배열 생성
heatmap_data = cell(length(y_range), length(x_range));

for r=1:height(all_data)
   x = all_data{r, 1};
   y = all_data{r, 2};
   value = (5 - all_data{r, 3}) / 4;
   
   x_idx = find(x_range == x);
   y_idx = find(y_range == y);
   
   % Append rt to the cell
   heatmap_data{y_idx, x_idx} = [heatmap_data{y_idx, x_idx}; value];
end

% Compute the mean of each cell
for i=1:numel(heatmap_data)
    heatmap_data{i} = mean(heatmap_data{i});
end

reduced_heatmap_data = zeros(length(y_range)/4, length(x_range)/4);
for i = 1:length(y_range)/4
    for j = 1:length(x_range)/4
        % Extract the 4x4 block
        block = cell2mat(heatmap_data((i-1)*4+1:i*4, (j-1)*4+1:j*4));
        % Replace 0 with NaN
        block(block == 0) = NaN;
        % Compute the mean ignoring NaN values
        reduced_heatmap_data(i, j) = nanmean(nanmean(block));
    end
end

imagesc(x_range(1:4:end), y_range(1:4:end), reduced_heatmap_data);
colormap('jet')
colorbar
xlabel('Potential Gain ($)');
ylabel('Potential Loss ($)');

% 제목 설정
title('Probability of Acceptance');
