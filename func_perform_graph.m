function graph = func_perform_graph(event_table, c_sbj)
    hold on
    colororder({'#0072BD','#000000'})

    % RT plot : left
    X=1:height(event_table);
    yyaxis left
    title([c_sbj ': RT & Correctness'],"FontSize",18,"FontWeight","bold")

    % RT plot
    plot(event_table.RT,'Color','#0072BD', 'LineWidth',1.8,'Marker','.','MarkerSize',20);
    % Timeout threshold
    yline(1.5,'-.','Timeout (>= 1.5s)','LabelHorizontalAlignment', 'center' ,'Color',"#0072BD",'LineWidth',1.2);

    xlabel('Trial','FontSize',15,'FontWeight','bold')
    ylabel('RT','FontSize',15,'FontWeight','bold')

    yticks(0:0.2:1.8)
    xlim([1 height(event_table)]);
    ylim([-0.2 2]);
    pbaspect([2 1 1]);

    % correctness plot : right
    yyaxis right
    Y=event_table.Correct_Num';
    Y(1,Y == 2) = 1;

    stairs(Y);
    coloringX = [X;X];
    coloringY = [Y;Y];
    area(coloringX([2:end end]),coloringY(1:end), 'FaceColor', 'k');
    ylim([0 10]);

    %legend
    img = imread('correctness_legend.jpg');
    image(img,'XData',[32 34],'YData',[3 0],'Clipping','off')

    %%%
    coloringTO = [event_table.isTimeout'; event_table.isTimeout'];
    area(coloringX([2:end end]), coloringTO(1:end), 'FaceColor',"#EDB120");
    yticks([])
    
    graph = gcf;
end
