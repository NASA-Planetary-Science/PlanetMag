function MakeHodograms

    figDir = 'figures';
    figXtn = 'pdf';
    figNumBase = 30;
    LIVE_PLOTS = 0;
    datDir = 'out';

    %% Jupiter family
    parentName = 'Jupiter';
    figNumber = figNumBase + 1;
    moons = [
        "Io"
        "Europa"
        "Ganymede"
        "Callisto"
        ];
    cfmt = {
        'r'
        'b'
        'g'
        'c'
        };
    nMoons = length(moons);

    windowName = [parentName ' system hodograms'];
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
    tiles = tiledlayout(1, nMoons, 'TileSpacing', 'tight', 'Padding', 'tight');

    for i=1:nMoons
        moonName = char(moons(i));
        if strcmp(moonName, 'Callisto')
            [~, ~, ~, ~, fEnd] = GetModelOpts(parentName, 6, -1);
            magModelDescrip = 'VIP4+K';
        else
            [~, ~, ~, ~, fEnd] = GetModelOpts(parentName, 0, -1);
            magModelDescrip = 'JRM33+C2020';
        end
        load(fullfile(datDir, ['evalB' moonName fEnd]), 'BvecMoon', 'outCoords')
        SPHOUT = 0;
        if ~strcmp(outCoords(1:3), 'IAU'); SPHOUT = 1; end

        nexttile
        PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, cfmt{i});
    end

    fName = [parentName 'FamilyHodogram'];
    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    fig.Units = fig.PaperUnits;
    fig.PaperSize = [1.025*fig.Position(3) 0.375*fig.Position(4)];
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end

    %% Saturn family
    parentName = 'Saturn';
    figNumber = figNumBase + 2;
    moons = [
        "Mimas"
        "Enceladus"
        "Dione"
        "Rhea"
        "Titan"
        ];
    cfmt = {
        'm'
        'c'
        'b'
        'k'
        '#f97306' % Orange
        };
    nMoons = length(moons);

    windowName = [parentName ' system hodograms'];
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
    ApplyPlotDefaults(fig, interpreter, font);
    tiles = tiledlayout(1, nMoons, 'TileSpacing', 'tight', 'Padding', 'tight');

    for i=1:nMoons
        moonName = char(moons(i));
        [~, ~, ~, ~, fEnd] = GetModelOpts(parentName, 0, -1);
        magModelDescrip = 'Cassini 11';
        load(fullfile(datDir, ['evalB' moonName fEnd]), 'BvecMoon', 'outCoords')
        SPHOUT = 0;
        if ~strcmp(outCoords(1:3), 'IAU'); SPHOUT = 1; end

        nexttile
        PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, cfmt{i});
    end

    fName = [parentName 'FamilyHodogram'];
    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    fig.Units = fig.PaperUnits;
    fig.PaperSize = [1.025*fig.Position(3) 0.3*fig.Position(4)];
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end

    %% Uranus family
    parentName = 'Uranus';
    figNumber = figNumBase + 3;
    moons = [
        "Miranda"
        "Ariel"
        "Umbriel"
        "Titania"
        "Oberon"
        ];
    cfmt = {
        'm'
        'r'
        'k'
        'g'
        '#06c2ac' % Turqoise
        };
    nMoons = length(moons);

    windowName = [parentName ' system hodograms'];
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
    ApplyPlotDefaults(fig, interpreter, font);
    tiles = tiledlayout(1, nMoons, 'TileSpacing', 'tight', 'Padding', 'tight');

    for i=1:nMoons
        moonName = char(moons(i));
        [~, ~, ~, magModelDescrip, fEnd] = GetModelOpts(parentName, 0, -1);
        load(fullfile(datDir, ['evalB' moonName fEnd]), 'BvecMoon', 'outCoords')
        SPHOUT = 0;
        if ~strcmp(outCoords(1:3), 'IAU'); SPHOUT = 1; end

        nexttile
        PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, cfmt{i});
    end

    fName = [parentName 'FamilyHodogram'];
    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    fig.Units = fig.PaperUnits;
    fig.PaperSize = [1.025*fig.Position(3) 0.3*fig.Position(4)];
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end

    %% Neptune family
    parentName = 'Neptune';
    figNumber = figNumBase + 4;
    cfmt = '#9a0eea'; % Violet


    moonName = 'Triton';
    [~, ~, ~, magModelDescrip, fEnd] = GetModelOpts(parentName, 0, -1);
    load(fullfile(datDir, ['evalB' moonName fEnd]), 'BvecMoon', 'outCoords')
    SPHOUT = 0;
    if ~strcmp(outCoords(1:3), 'IAU'); SPHOUT = 1; end
    titleInfo = [moonName ', ' magModelDescrip];

    fName = [parentName 'FamilyHodogram'];
    windowName = [parentName ' system hodograms'];
    [xx, yy, ~, ~, xInfo, yInfo, ~, xlims, ylims] = HodogramSingle( ...
        moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, 0, 0);
    PlotGeneric(xx, yy, [], windowName, titleInfo, xInfo, yInfo, fName, figDir, figXtn, ...
        LIVE_PLOTS, figNumber, 'linear', 'linear', xlims, ylims, cfmt, 1, 1);
end

    function [xInfo, yInfo] = PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, ...
            SPHOUT, cfmt)
    [xx, yy, ~, ~, xInfo, yInfo, ~, xlims, ylims] = HodogramSingle( ...
            moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, 0, 0);

    plot(xx, yy, 'color', cfmt);
    pbaspect([1 1 1])
    title([moonName ', ' magModelDescrip]);
    xlim(xlims);
    ylim(ylims);
end