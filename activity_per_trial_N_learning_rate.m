%% activity and learning rate DTW
close all
% load('regressors_GLM_0922.mat','roi_hpc_name','roi_ctx_name')
cd('G:\JSR\240922_new_fmri\Activity')
load('G:\JSR\240922_new_fmri\bhv\data_learning_curve\curve_95/pdata_table.mat')
load('../learn_curv_info.mat')
acq_95_41(exclud_sbj)=[];

dtw_output = struct;
max_len = 0;

hpc_or_ctx={'ctx'};     %{'hpc','ctx'}

curr_names=fieldnames(activ.(hpc_or_ctx{:}));
curr_names=curr_names(~strcmp(curr_names,'trials'));

aligned_info = table;

%%
for perform_group={'all','early','late','fail'}
plot_count = 0;
figure_index = 1;

if strcmp(perform_group,'all')
    idx=1:numel(sbj_id_list_41);    
elseif strcmp(perform_group,'early')
    idx=idx_early;    
elseif strcmp(perform_group,'late')
    idx=idx_late;    
elseif strcmp(perform_group,'fail')
    idx=idx_fail;
end
    sbj_group=sbj_id_list_41(idx);
    sbj_learn_trial=acq_95_41(idx);

for r=1:numel(curr_names)
    
    % First pass to determine the maximum length
    for sbj_i = 1:numel(sbj_group)
        y = pdata_table.(sprintf('sub_%.2d', sbj_group(sbj_i)))';
        x = activ.(hpc_or_ctx{:}).(curr_names{r})(idx(sbj_i), :);
        
        [~, ix, iy] = dtw(x, y);
        max_len = max(max_len, length(ix));
        max_len = max(max_len, length(iy));
    end
    
    aligned_x_all = nan(numel(sbj_group), max_len);
    aligned_y_all = nan(numel(sbj_group), max_len);
    aligned_learn_trial = nan(numel(sbj_group), 1);
    y_mat=[];
    p_values = nan(numel(sbj_group), 1); 
    
    for sbj_i = 1:numel(sbj_group)
        y = pdata_table.(sprintf('sub_%.2d',sbj_group(sbj_i)))';
        y_mat(sbj_i,:)=y;
        x=activ.(hpc_or_ctx{:}).(curr_names{r})(idx(sbj_i),:);
        
        [dist, ix, iy] = dtw(x, y);
        
        diff_t = x(ix) - y(iy);
        
        [h, p] = ttest(diff_t);
        p_values(sbj_i) = p; % p 값을 배열에 저장
        
        dtw_output.(sprintf('sub_%.2d', sbj_group(sbj_i))) = {'dist', dist; 'ix', ix'; 'iy', iy'; 'p-value', p};
        
        aligned_x_all(sbj_i, 1:length(ix)) = x(ix);
        aligned_y_all(sbj_i, 1:length(iy)) = y(iy);
        
        if ~isempty(find(iy==sbj_learn_trial(sbj_i),1))
        aligned_learn_trial(sbj_i) = find(iy==sbj_learn_trial(sbj_i),1); 
        end
    end
    
    mean_aligned_x = nanmean(aligned_x_all, 1);
    mean_aligned_y = nanmean(aligned_y_all, 1);
    mean_p_value = nanmean(p_values); 
    if ~isempty(aligned_learn_trial)
        mean_trial = nanmean(aligned_learn_trial, 1); 
    end
    
    %% plotting
    %% plotting original signals
    if mod(plot_count, 6) == 0
        figure('position',[827,536,1390,725]);
        figure_index = figure_index + 1;
    end
    
    plot_count = plot_count + 1;
    
    % Plot original signals
    subplot(3, 2, mod(plot_count-1, 6) + 1);
    yyaxis left;
    plot(mean(activ.(hpc_or_ctx{:}).(curr_names{r}), 1), '-o');
    ylabel('activity');
    
    yyaxis right;
    plot(mean(y_mat, 1), '-x');
    ylabel('learning rate');
    title(sprintf('Mean Original Signals (%s - %s)', hpc_or_ctx{:}, strrep(curr_names{r}, '_', '.')));
    
    if mod(plot_count, 3) == 0
        legend('activity', 'learning rate','Position',[0.889172287757767,0.015557516378381,0.094834884150769,0.050065018188194]);
    elseif mod(plot_count, 6) == 4
        legend('activity', 'learning rate','Position',[0.887014014376472,0.611419585343898,0.094834884150769,0.050065018188194]);
    end
    
    %% plotting aligned signals
    if mod(plot_count, 6) == 0
        figure('position',[827,536,1390,725]);
        figure_index = figure_index + 1;
    end
    
    plot_count = plot_count + 1;
    
    subplot(3, 2, mod(plot_count-1, 6) + 1);
    yyaxis left;
    if mean_p_value < 0.1
        plot(mean_aligned_x, '-o', 'LineWidth', 2); % mean p value가 0.1 이하일 때 선 굵게 표시
    else
        plot(mean_aligned_x, '-o');
    end
    ylabel('activity (aligned)');
    
    yyaxis right;
    if mean_p_value < 0.1
        plot(mean_aligned_y, '-x', 'LineWidth', 2); % mean p value가 0.1 이하일 때 선 굵게 표시
    else
        plot(mean_aligned_y, '-x');
    end
    ylabel('learning rate (aligned)');
    title(sprintf('Mean Aligned Signals (%s - %s)', hpc_or_ctx{:}, strrep(curr_names{r}, '_', '.')));
    
    % mean_learning trial을 세로로 그리기
if exist('mean_trial', 'var') && ~isnan(mean_trial) && ~isempty(mean_trial)
    xline(mean_trial, '--k', 'LineWidth', 1.5); 
end
    
     % mean_p_value < 0.1일 때 align 정보 저장
    if mean_p_value < 0.1
        for k = 1:length(mean_ix)
            new_row = {hpc_or_ctx{:}, curr_names{r}, perform_group{:}, mean_ix(k), mean_iy(k)};
            aligned_info = [aligned_info; new_row];
        end
    end
    
    save_folder=sprintf('G:/JSR/240922_new_fmri/Activity/DTW/%s_sub',perform_group{:});
    mkdir(save_folder)
    saveas(gcf, sprintf('%s/DTW_%s_%d.png', save_folder,hpc_or_ctx{:},figure_index-1));
    
end
end

aligned_info.Properties.VariableNames = {'hpc_cortex', 'ROI', 'perform_group', 'trial_neural_activity', 'trial_learn'};

aligned_info

writetable(aligned_info,'./DTW/aligned_info.xlsx')