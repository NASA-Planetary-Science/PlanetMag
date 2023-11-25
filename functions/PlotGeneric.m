function PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, figDir, fName, ...
    figXtn, LIVE_PLOTS)
    
    if LIVE_PLOTS
        figVis = 'on';
    else
        figVis = 'off';
    end

    fig = figure('Visible', figVis); hold on;
    set(gcf,'Name', windowName);
    yySize = size(yy);
    nPlots = yySize(1);
    for i=1:nPlots
        plot(xx, yy(i,:));
        plot(xx, yy(i,:));
    end
    title(titleInfo);
    xlabel(xInfo);
    ylabel(yInfo);
    legend(legendStrings);
    if ~LIVE_PLOTS
        outFig = fullfile(figDir, [fName '.' figXtn]);
        saveas(fig, outFig)
        disp(['Figure saved to ' outFig '.'])
    end
end