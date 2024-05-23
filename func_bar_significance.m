function func_bar_significance(b, p, err)

for b_i=1:max(size(b))
    curr_ax = b(b_i).Parent;
    curr_b=b(b_i);
    curr_err=err(b_i);
    curr_p=p(b_i);


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

    hold off

end
end