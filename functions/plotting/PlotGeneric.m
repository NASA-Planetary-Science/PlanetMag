function fig = PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, ...
    fName, figDir, figXtn, LIVE_PLOTS, figNumber, xAxisScale, yAxisScale, xlims, ylims, cfmt, ...
    DO_BOX, DO_GRID, vlinexCols, ticksx, ticklabelsx, ticksy, ticklabelsy)
% Make a simple plot of data and optionally print to disk.
%
% Parameters
% ----------
% xx : double, 1xN
%   x axis data.
% yy : double, MxN
%   y axis data. Each row will be plotted and must have the same number of columns as the length of
%   ``xx``.
% legendStrings : string, 1xM'
%   A list of strings with length matching the number of rows in ``yy`` (or fewer) with which to
%   describe the plotted data in the legend. If the number of strings passed here is fewer than the
%   number of rows in ``yy``, only the first M' plots will be labeled in the legend.
% windowName : char, 1xC
%   Name to use for the interactive figure window, which is shown when ``LIVE_PLOTS`` is true.
% titleInfo : char, 1xD
%   Description of plot contents to print at the top of the figure.
% xInfo : char, 1xE
%   Label for x axis of plot.
% yInfo : char, 1xF
%   Label for y axis of plot.
% fName : char, 1xG
%   File name to use for output figures, sans extension.
% figDir : char, 1xH, default='figures'
%   Directory to use for output figures.
% figXtn : char, 1xI, default='pdf'
%   Extension to use for figures, which determines the file type.
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% figNumber : int, default=0
%   Figure number to use for assigning/reusing figure windows. If ``0`` is passed, the default (no
%   assigned number, use increment) is used.
% xAxisScale : char, default='linear'
%   Scale type of the x axis, either ``'linear'`` or ``'log'``.
% yAxisScale : char, default='linear'
%   Scale type of the y axis, either ``'linear'`` or ``'log'``.
% xlims : double, 1x2, default='auto'
%   Minimum and maximum values to display on x axis, respectively, or the keyword ``'auto'``.
% ylims : double, 1x2, default='auto'
%   Minimum and maximum values to display on y axis, respectively, or the keyword ``'auto'``.
% cfmt : cell, 1xM, default={}
%   Color/style format descriptions. Must contain a cell for each row in ``yy``. Cells may contain
%   a color/style format string (e.g. '--k' for a dashed black line) or other color-setting options
%   such as an RGB vector.
% DO_BOX : bool, default=0
%   Add an outer box to plots.
% DO_GRID : bool, default=0
%   Add an outer box to plots.
% vlinexCols : cell, default={}
%   Add vertical lines with given color to plots. If passed, vlinexCols{1} must be a vector of x
%   positions in data units and vlinexCols{2} must be a color.
% ticksx : double, 1xU, default=[]
%   Explicit values at which to place tick marks on the x axis.
% ticklabelsx : cell, 1xU, default={}
%   Char array labels for tick marks on the x axis passed in ``ticksx``.
% ticksy : double, 1xV, default=[]
%   Explicit values at which to place tick marks on the y axis.
% ticklabelsy : cell, 1xV, default={}
%   Char array labels for tick marks on the y axis passed in ``ticksy``.
%
% Returns
% -------
% fig : figure
%   Figure object for the generated plot.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('figNumber', 'var'); figNumber = 0; end
    if ~exist('xAxisScale', 'var'); xAxisScale = 'linear'; end
    if ~exist('yAxisScale', 'var'); yAxisScale = 'linear'; end
    if ~exist('xlims', 'var'); xlims = 'auto'; end
    if ~exist('ylims', 'var'); ylims = 'auto'; end
    if ~exist('cfmt', 'var'); cfmt = {}; end
    if ~exist('DO_BOX', 'var'); DO_BOX = 0; end
    if ~exist('DO_GRID', 'var'); DO_GRID = 0; end
    if ~exist('vlinexCols', 'var'); vlinexCols = {}; end
    if ~exist('ticksx', 'var'); ticksx = []; end
    if ~exist('ticklabelsx', 'var'); ticklabelsx = {}; end
    if ~exist('ticksy', 'var'); ticksy = []; end
    if ~exist('ticklabelsy', 'var'); ticklabelsy = {}; end

    if LIVE_PLOTS
        if ~figNumber
            fig = figure('Visible', 'on', 'Name', windowName);
        else
            fig = figure(figNumber);
            clf();
            set(gcf, 'Visible', 'on', 'Name', windowName);
        end
    else
        fig = figure('Visible', 'off', 'Name', windowName);
    end
    hold on;
    [interpreter, font] = SetPlotDefaults();
    ApplyPlotDefaults(fig, interpreter, font);

    yySize = size(yy);
    nPlots = yySize(1);
    DO_COLORS = length(cfmt) == nPlots;
    if length(cfmt) > 1 && ~DO_COLORS
        warning(['Only ' num2str(length(cfmt)) ' colors for ' num2str(nPlots) ' were passed' ...
            'to PlotGeneric.'])
    end

    for i=1:nPlots
        if DO_COLORS
            plot(xx, yy(i,:), cfmt{i});
        else
            plot(xx, yy(i,:));
        end
    end
    if ~isempty(vlinexCols)
        vlines(vlinexCols{1}, vlinexCols{2});
    end
    if DO_BOX; box on; end
    if DO_GRID; grid on; end

    if ~isempty(ticksx); xticks(ticksx); end
    if ~isempty(ticklabelsx); xticklabels(ticklabelsx); end
    if ~isempty(ticksy); xticks(ticksy); end
    if ~isempty(ticklabelsy); yticklabels(ticklabelsy); end

    title(titleInfo);
    xlabel(xInfo);
    ylabel(yInfo);
    xlim(xlims);
    ylim(ylims);
    % Make axes square if x and y limits have same numbers
    if isnumeric(xlims) && isnumeric(ylims) && all(xlims == ylims)
        pbaspect([1 1 1])
    end
    set(gca, 'xscale', xAxisScale);
    set(gca, 'yscale', yAxisScale);
    if ~isempty(legendStrings); legend(legendStrings); end

    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    fig.Units = fig.PaperUnits;
    fig.PaperSize = fig.Position(3:4);
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end
end