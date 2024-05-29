
function func_bar_significance(b, p, err)

%% input
%    b= bar plot. but b doesn't have to be a 'bar' plot. 
%    p= output of ttest
%    err= output of errorbar

 %%
% for b_i=1:max(size(b))
%     curr_ax = b(b_i).Parent;
%     curr_b=b(b_i);
%     curr_err=err(b_i);
%     curr_p=p(b_i);
% 
%     hold(curr_ax ,'on');
% 
%     dy = 0.1 * curr_ax.YLim(2);
%     y_star = 0.8 * dy;
% 
%     idx_plus = curr_p < 0.1 & curr_p >= 0.05;
%     idx_star = curr_p < 0.05 & curr_p >= 0.01;
%     idx_double_star = curr_p < 0.01 & curr_p >= 0.001;
%     idx_triple_star = curr_p < 0.001;
% 
%     add_star = @(i, symbol) text(curr_ax, curr_b.XData(i), curr_err.YPositiveDelta(i) + curr_b.YData(i) + y_star, symbol, ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20, 'LineStyle', '-');
% 
%     if any(idx_plus)
%         for i = find(idx_plus)'
%             add_star(i, '+');
%         end
%     end
% 
%     if any(idx_star)
%         for i = find(idx_star)'
%             add_star(i, '*');
%         end
%     end
% 
%     if any(idx_double_star)
%         for i = find(idx_double_star)'
%             add_star(i, '**');
%         end
%     end
% 
%     if any(idx_triple_star)
%         for i = find(idx_triple_star)'
%             add_star(i, '***');
%         end
%     end
% 
%     hold off
% 
% end

%%
if numel(b)>1
curr_ax = b(1).Parent;
    curr_b=b(1);
    curr_err=err(1);
    curr_p=p(1,:);


    hold(curr_ax ,'on');

    dy = 0.1 * curr_ax.YLim(2);
    y_star = 0.8 * dy;

    idx_plus = curr_p < 0.1 & curr_p >= 0.05;
    idx_star = curr_p < 0.05 & curr_p >= 0.01;
    idx_double_star = curr_p < 0.01 & curr_p >= 0.001;
    idx_triple_star = curr_p < 0.001;

    add_star = @(i, symbol) text(curr_ax, curr_b.XData(i)-0.13, curr_err.YPositiveDelta(i) + curr_b.YData(i) + y_star, symbol, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20, 'LineStyle', '-','Color','#D95319');

    if any(idx_plus)
        for i = find(idx_plus)'
            add_star(i, '+');
        end
    end

    if any(idx_star)
        for i = find(idx_star)'
            add_star(i, '*');
        end
    end

    if any(idx_double_star)
        for i = find(idx_double_star)'
            add_star(i, '**');
        end
    end

    if any(idx_triple_star)
        for i = find(idx_triple_star)'
            add_star(i, '***');
        end
    end

%%

curr_ax = b(2).Parent;
    curr_b=b(2);
    curr_err=err(2);
    curr_p=p(2,:);


    hold(curr_ax ,'on');

    dy = 0.1 * curr_ax.YLim(2);
    y_star = 0.8 * dy;

    idx_plus = curr_p < 0.1 & curr_p >= 0.05;
    idx_star = curr_p < 0.05 & curr_p >= 0.01;
    idx_double_star = curr_p < 0.01 & curr_p >= 0.001;
    idx_triple_star = curr_p < 0.001;

    add_star = @(i, symbol) text(curr_ax, curr_b.XData(i)+0.13, curr_err.YPositiveDelta(i) + curr_b.YData(i) + y_star, symbol, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20, 'LineStyle', '-','Color','#0072BD');

    if any(idx_plus)
        for i = find(idx_plus)'
            add_star(i, '+');
        end
    end

    if any(idx_star)
        for i = find(idx_star)'
            add_star(i, '*');
        end
    end

    if any(idx_double_star)
        for i = find(idx_double_star)'
            add_star(i, '**');
        end
    end

    if any(idx_triple_star)
        for i = find(idx_triple_star)'
            add_star(i, '***');
        end
    end





    hold off


elseif numel(b)==1
    curr_ax = b.Parent;
    curr_b=b;
    curr_err=err;
    curr_p=p;


    hold(curr_ax ,'on');

    dy = 0.1 * curr_ax.YLim(2);
    y_star = 0.8 * dy;

    idx_plus = curr_p < 0.1 & curr_p >= 0.05;
    idx_star = curr_p < 0.05 & curr_p >= 0.01;
    idx_double_star = curr_p < 0.01 & curr_p >= 0.001;
    idx_triple_star = curr_p < 0.001;

    add_star = @(i, symbol) text(curr_ax, curr_b.XData(i), curr_err.YPositiveDelta(i) + curr_b.YData(i) + y_star, symbol, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20, 'LineStyle', '-');

    if any(idx_plus)
        for i = find(idx_plus)'
            add_star(i, '+');
        end
    end

    if any(idx_star)
        for i = find(idx_star)'
            add_star(i, '*');
        end
    end

    if any(idx_double_star)
        for i = find(idx_double_star)'
            add_star(i, '**');
        end
    end

    if any(idx_triple_star)
        for i = find(idx_triple_star)'
            add_star(i, '***');
        end
    end




end