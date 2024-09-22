data_dir='G:\JSR\spm_prep_glm_main_0804\final2\across_within_comparison';
cd(data_dir)
curr_names={'Lt_Hp','Lt_DGCA3','Lt_CA1','Rt_Hp','Rt_DGCA3','Rt_CA1','Bi_Hp','Bi_DGCA3','Bi_CA1'};

%% roi만 저장
same_pat_xls='main_ver_2_within_same_hpc_lap_final_pattern_0827.xlsx';
diff_pat_xls='main_ver_2_within_diff_hpc_lap_final_pattern_0827.xlsx';


for roi_i=1:numel(curr_names)
    curr_roi=curr_names{roi_i};
    same=readtable(same_pat_xls,'Sheet',curr_roi);
    diff=readtable(diff_pat_xls,'Sheet',curr_roi);
    
    writetable(same, 'same_roi_lap_pattern_0905.xlsx', 'Sheet',curr_roi);
    writetable(diff, 'diff_roi_lap_pattern_0905.xlsx', 'Sheet',curr_roi);

end

%% struct 제작 from ROI pattern per lap

same_pat_xls='same_roi_lap_pattern_0905.xlsx';
diff_pat_xls='diff_roi_lap_pattern_0905.xlsx';

roi_pat = struct;
for roi_i=1:numel(curr_names)
    curr_roi=curr_names{roi_i};
    
    same=readtable(same_pat_xls,'Sheet',curr_roi);
    diff=readtable(diff_pat_xls,'Sheet',curr_roi);
    
        roi_pat.same.(curr_roi)=table2array(same);
        roi_pat.diff.(curr_roi)=table2array(diff);

%     for l=1:8
%         roi_pat.same.(curr_roi).(sprintf('lap_%d',l))=table2array(same(:,l));
%         roi_pat.diff.(curr_roi).(sprintf('lap_%d',l))=table2array(diff(:,l));
%     end
end
        
  save('roi_pat_laps.mat',"roi_pat")
        
        
%% %%%%%%%%%% 240913 %%%%%%%%%%%%%%%%%
%% learning rate 구하기
parray_34=[];
for sbj_i=1:34
    c_sbj=sprintf('sub_%.2d',sbj_id_list_38(sbj_i));
    
    parray_34(:,sbj_i)=ptable.(c_sbj);
end
% parray_34=parray_34';
%% learning rate per lap 
parray_lap=[];
for l=1:8
parray_lap(:,l) = mean(parray_34((l-1)*4+1:l*4, :), 1);
end

%% correlation btwn learning rate and roi pattern 

for roi_i=1:numel(curr_names)
    curr_roi=curr_names{roi_i};
    s=roi_pat.same.(curr_roi);
    d=roi_pat.diff.(curr_roi);
    
    for l = 1:8
        [corr_lap.same(l), p_values.same(l)] = corr(s(:, l), parray_lap(:, l));
        [corr_lap.diff(l), p_values.diff(l)] = corr(d(:, l), parray_lap(:, l));
    end
end

%same
corrplot([s, parray_lap])
% find(temp_p_values<0.1)

[corr_same,p_v_same]=corr([s, parray_lap]);

% idx=find(p_v_same<0.1);
% sig_p=zeros(size(p_v_same));
% sig_p(idx)=p_v_same(idx);
% sig_corr=zeros(size(p_v_same));
% sig_corr(idx)=corr_same(idx);

sig_p = p_v_same;
sig_corr = corr_same;

sig_p(p_v_same >= 0.1) = NaN;
sig_corr(p_v_same >= 0.1) = NaN;

%diff
corrplot([d, parray_lap])
% find(temp_p_values<0.1)

[corr_diff,p_v_diff]=corr([d, parray_lap]);

sig_p = p_v_diff;
sig_corr = corr_diff;

% Set values with p >= 0.1 to zero
sig_p(p_v_diff >= 0.1) = NaN;
sig_corr(p_v_diff >= 0.1) = NaN;





% 
% for i = 1:8
%     for j = 1:8
%         [temp_cor(i, j), temp_p_values(i, j)] = corr(s(:, i), parray_lap(:, j));
%     end
% end

%% p-value and plot        

signif_lap_same=find(p_values.same<0.1);
signif_lap_diff=find(p_values.diff<0.1);

corr(s(:, signif_lap_same), parray_lap(:, signif_lap_same))
corr(d(:, signif_lap_diff), parray_lap(:, signif_lap_diff))





figure;
subplot(2, 1, 1);
bar(corr_lap.diff);
title('Correlation Coefficients');
xlabel('Column Index');
ylabel('Correlation Coefficient');

subplot(2, 1, 2);
bar(p_values.diff);
title('P-values');
xlabel('Column Index');
ylabel('P-value');
        
        