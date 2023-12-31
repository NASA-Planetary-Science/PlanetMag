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
        '#bf9005' % ochre
        '#0652ff' % electric blue
        '#76cd26' % apple green
        '#75bbfd' % sky blue
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
        '#c20078' % magenta
        '#8e82fe' % periwinkle
        '#02ab2e' % kelly green
        '#380835' % eggplant
        '#dbb40c' % gold
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
        PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, cfmt{i}, 0);
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
        '#650021' % maroon
        '#13eac9' % aqua
        'k' % black
        '#06c2ac' % turqoise
        '#f97306' % orange
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
    cfmt = '#d1b26f'; % tan

    moonName = 'Triton';
    [~, ~, ~, magModelDescrip, fEnd] = GetModelOpts(parentName, 0, -1);
    load(fullfile(datDir, ['evalB' moonName fEnd]), 'BvecMoon', 'outCoords')
    SPHOUT = 0;
    if ~strcmp(outCoords(1:3), 'IAU'); SPHOUT = 1; end
    titleInfo = [moonName ', ' magModelDescrip];

    fName = [parentName 'FamilyHodogram'];
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
    nexttile
    PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, cfmt);

    fName = [parentName 'FamilyHodogram'];
    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    fig.Units = fig.PaperUnits;
    fig.PaperSize = [0.25*fig.Position(3) 0.35*fig.Position(4)];
    fig.PaperPosition(1:2) = [0.05 -1.5];
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end
end

    function [xInfo, yInfo] = PlotHodogramTile(moonName, parentName, magModelDescrip, BvecMoon, ...
            SPHOUT, cfmt, LIMS)
    if ~exist('LIMS', 'var'); LIMS = 1; end
    [xx, yy, ~, ~, xInfo, yInfo, ~, xlims, ylims] = HodogramSingle( ...
            moonName, parentName, magModelDescrip, BvecMoon, SPHOUT, 0, 0);

    plot(xx, yy, 'color', cfmt);
    pbaspect([1 1 1])
    daspect([1 1 1])
    title([moonName ', ' magModelDescrip]);
    if LIMS
        xlim(xlims);
        ylim(ylims);
    end
end