function MakeHodograms
% Print hodograms for each moon implemented in PlanetProfile from the default model data used to
% calculate the excitation moments.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figDir = 'figures';
    figXtn = 'pdf';
    figNumBase = 30;
    LIVE_PLOTS = 0;
    datDir = 'out';
    tileSize = 2; % in inches
    pad = 0.025*tileSize;

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
        '#ec2d01' % tomato red
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
    figWidth = nMoons*tileSize+2*pad;
    figHeight = tileSize+2*pad;
    fig.Units = 'inches';
    fig.Position(3:4) = [nMoons*tileSize figHeight];
    fig.PaperSize = [figWidth figHeight];
    fig.PaperPositionMode = 'auto';
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end

    %% Saturn family
    parentName = 'Saturn';
    % Part 1
    figNumber = figNumBase + 2;
    moons = [
        "Mimas"
        "Enceladus"
        "Tethys"
        "Dione"
        ];
    cfmt = {
        '#c20078' % magenta
        '#8e82fe' % periwinkle
        '#fe420f' % orangered
        '#02ab2e' % kelly green
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

    fName = [parentName 'FamilyHodogram1'];
    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    figWidth = nMoons*tileSize+2*pad;
    fig.Units = 'inches';
    fig.Position(3:4) = [nMoons*tileSize figHeight];
    fig.PaperSize = [figWidth figHeight];
    fig.PaperPositionMode = 'auto';
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end

    % Saturn part 2
    figNumber = figNumBase + 5;
    moons = [
        "Rhea"
        "Titan"
        "Iapetus"
        ];
    cfmt = {
        '#9900fa' % vivid purple
        '#dbb40c' % gold
        '#fc5a50' % coral
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

    fName = [parentName 'FamilyHodogram2'];
    outFig = fullfile(figDir, [fName '.' figXtn]);
    % Crop page size for pdf outputs
    figWidth = nMoons*tileSize+2*pad;
    fig.Units = 'inches';
    fig.Position(3:4) = [nMoons*tileSize figHeight];
    fig.PaperSize = [figWidth figHeight];
    fig.PaperPositionMode = 'auto';
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
        '#8c000f' % crimson
        '#5cac2d' % grass
        'k' % black
        '#be03fd' % bright purple
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
    figWidth = nMoons*tileSize+2*pad;
    fig.Units = 'inches';
    fig.Position(3:4) = [nMoons*tileSize figHeight];
    fig.PaperSize = [figWidth figHeight];
    fig.PaperPositionMode = 'auto';
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end

    %% Neptune family
    parentName = 'Neptune';
    figNumber = figNumBase + 4;
    cfmt = '#78d1b6'; % seafoam blue
    nMoons = 1;

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
    figWidth = nMoons*tileSize+2*pad;
    fig.Units = 'inches';
    fig.Position(3:4) = [nMoons*tileSize figHeight];
    fig.PaperSize = [figWidth figHeight];
    fig.PaperPositionMode = 'auto';
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