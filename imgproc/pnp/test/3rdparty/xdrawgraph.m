function xdrawgraph(xs,yrange,method_list,field,ti,lx,ly)
    %the legend is at upper right in default
    box('on');
    hold('all');

    p= zeros(size(method_list));
    lgs = [];
    for i= 1:length(method_list)
        p(i)= plot(xs,method_list(i).(field),'marker',method_list(i).marker,...
            'color',method_list(i).color,...
            'markerfacecolor',method_list(i).markerfacecolor,...
            'displayname',method_list(i).name, ...
            'LineWidth',2,'MarkerSize',8,'LineStyle','-');
        lgs = [lgs, method_list(i).name];
    end
    ylim(yrange);
    xlim([0 xs(end)]);
    set(gca,'xtick',xs);

    title(ti,'FontSize',12,'FontName','Arial','Interpreter','latex');
    xlabel(lx,'FontSize',11,'Interpreter','latex');
    ylabel(ly,'FontSize',11,'Interpreter','latex');
    legend(p,lgs,'Interpreter','latex');
    grid on;
end