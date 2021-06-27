function PlotConcludeFigures( SET , OPTIONS , nameOfOptFunction, MinCost  )
% Plot some results
% link to the num of OptFunction and Bench in Monte in order to plot sub graph.

m = SET.row_of_subplot ;
n = SET.col_of_subplot ;
k = SET.index_of_subplot;


switch nameOfOptFunction
    case 'BBO'
        subplot(m , n , k);
        plot([0:OPTIONS.Maxgen], MinCost, '-k')
        title('BBO')
        xlabel('Generation');
        ylabel('Minimum Cost');
    case 'BBO_PI'
        subplot(m , n , k);
        plot([0:OPTIONS.Maxgen], MinCost, ':r')
        title('BBO by PI')
        xlabel('Generation');
        ylabel('Minimum Cost');
    case 'BBO_TI'
        subplot(m , n , k);
        plot([0:OPTIONS.Maxgen], MinCost, '-r')
        title('BBO by TI')
        xlabel('Generation');
        ylabel('Minimum Cost');
    case 'BBO_PE'
        subplot(m , n , k);
        plot([0:OPTIONS.Maxgen], MinCost, ':b')
        title('BBO by PE')
        xlabel('Generation');
        ylabel('Minimum Cost');
    case 'BBO_TE'
        subplot(m , n , k);
        plot([0:OPTIONS.Maxgen], MinCost, '-b')
        title('BBO by TE')
        xlabel('Generation');
        ylabel('Minimum Cost');
end

end

