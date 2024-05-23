clear;clc;close all
code_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code'; addpath(code_path)
path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\new_ocat_glm\ZY\figures\240523';
addpath(genpath(path)); cd(path)
load('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code\odt_data.mat');

sbj_list(sbj_list==4)=[]; %4번은 ODT parsing이 잘못되어있었음!
fail_sbj=[5 6 15 22 30];

%% input

which_sbj='learn'; %'all', 'learn', 'fail'
roi = {'Lt.Hp', 'Lt.PHC', 'Lt.CA23DG', 'Lt.CA1', ...
    'Rt.Hp', 'Rt.PHC',  'Rt.CA23DG','Rt.CA1', ...
    'Bi.Hp', 'Bi.PHC',  'Bi.CA23DG','Bi.CA1'};

% sbj indexing
if strcmp(which_sbj,'learn')
    sbj_idx= find(~ismember(sbj_list, [fail_sbj 4]));
    ttl = "learned subjects (n=22)";

elseif strcmp(which_sbj,'fail')
    sbj_idx= find(ismember(sbj_list, [fail_sbj 4]));
    ttl = "failed subjects (n=5)";
else
    sbj_idx= find(sbj_list~=4);
    ttl = "All subjects(n=27)";
end

% region indexing
region_idx=find(ismember(hpc_name_all, roi));
hpc_region=hpc_name_all(region_idx);


same=[];diff=[];
for sbj_i=sbj_idx
    same=[same; same_pre_post_diff{1, sbj_i}(region_idx)];
    diff=[diff; diff_pre_post_diff{1, sbj_i}(region_idx)];
end

writematrix(same, "ps_change.xlsx", 'Sheet', 'same');
writematrix(diff, "ps_change.xlsx", 'Sheet', 'diff');

%% samd and diff context
M=[mean(same); mean(diff)]';

figure('Position',[246,361,1498,952]);
ax = gca;
h = bar(ax,M, 'grouped','LineWidth',1,'FaceColor','flat');
hold(ax,'on');
err=errorbar(ax,cat(1,h(:).XEndPoints)',cat(1,h(:).YEndPoints)',0.05.*rand(size(h(1).YData)),'Color','k','LineWidth',1.5,'LineStyle','none');

set(ax, 'XTick', 1:size(M,1), 'XTickLabel', hpc_name_all(region_idx),'LineWidth',1, 'FontSize',18,'FontWeight','bold');
xtickangle(45); h(1).FaceColor = '#D95319'; h(2).FaceColor = '#0072BD';legend(h([1 2]), {'same', 'diff'}); box off; title(ttl)

%stats
[~, p_same] = ttest(same);
[~, p_diff] = ttest(diff);
func_bar_significance(h, [p_same,p_diff] , err)

hold on;
[~, p] = ttest2(same, diff);
func_bar_significance(h(1), p, err)



%% samd and diff context_individual

figure('position', [-1070,55,3570,1457])

% for문 시작
for i=1:height(same)
    % subplot 추가
    subplot(4,6,i)
    h=bar([same(i,:); diff(i,:)]','grouped','LineWidth',1,'FaceColor','flat')

    ax = gca;  hold(ax,'on');
    set(ax, 'XTick', 1:length(same), 'XTickLabel', hpc_name_all(region_idx),'LineWidth',1, 'FontSize',12,'FontWeight','bold', 'YLim', [-0.5 0.5],'YTick', -0.5:0.1:0.5);
    xtickangle(45); h(1).FaceColor = '#D95319'; h(2).FaceColor = '#0072BD'; box off; title(sprintf('subject# %.2d',sbj_list(sbj_idx(i))))

end
hold(ax,'off')
hold on
legend(h([1 2]), {'same', 'diff'},'Location','south')
hold off

saveas(gcf,strcat(which_sbj,'_indiv_Pattern_similarity_change.jpg'))


%% same-diff difference_all
stbr=same-diff;
[~, p] = ttest(stbr);
writematrix(stbr, "ps_change.xlsx", 'Sheet', 'same-diff');

figure('Position',[246,361,1498,952]);
b = bar(mean(stbr), 'LineWidth', 1,'FaceColor',[0.8 0.8 0.8]);
ax = gca; hold(ax,'on');
err = errorbar(ax, b.XEndPoints, b.YEndPoints, 0.05.*rand(size(b.YData)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', 'none');
set(ax, 'XTick', 1:numel(region_idx), 'XTickLabel', hpc_name_all(region_idx), 'LineWidth', 1, 'FontSize', 18, 'FontWeight', 'bold');
xtickangle(45); title(ttl);box off;

func_bar_significance(b, p,err)


saveas(ax,strcat(which_sbj,'_','Pattern_same-diff.jpg'))



%% same-diff difference_individual

figure('position', [-1070,55,3570,1457])

for i=1:height(stbr)
    % subplot 추가
    subplot(4,6,i)

 b=bar(stbr(i,:),'LineWidth',1,'FaceColor',[0.8 0.8 0.8])

    ax = gca;  hold(ax,'on');
    set(gca, 'XTick', 1:numel(region_idx), 'XTickLabel', hpc_name_all(region_idx),'LineWidth',1, 'FontSize',12,'FontWeight','bold', 'YLim', [-0.5 0.5],'YTick', -0.5:0.1:0.5);
    xtickangle(45); box off; title(sprintf('subject# %.2d',sbj_list(sbj_idx(i))))
    hold(ax,'off')
end
saveas(gcf,strcat(which_sbj,'_indiv_Pattern_same-diff.jpg'))












