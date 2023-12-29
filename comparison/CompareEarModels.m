function CompareEarModels(LIVE_PLOTS, scName, xtn, SEQUENTIAL, coeffPath, figDir, figXtn)
% Compare magnetic field measurements from spacecraft near Earth against each implemented magnetic
% field model.
%
% Datasets:
%
%   * Swarm MAG data: https://earth.esa.int/web/guest/swarm/data-access
%       * Current direct link for 1 s decimated measurements:
%         https://swarm-diss.eo.esa.int/#swarm/Level1b/Latest_baselines/MAGx_LR
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% scName : string, default="Swarm"
%   Spacecraft name for which magnetic field data will be compared against implemented models. A
%   directory must exist with this name in the ``MAG`` directory within the top-level P\lanetMag
%   directory. This directory will be searched for data files with the ``xtn`` extension, and each
%   of these files will be loaded.
% xtn : char, 1xC, default='.tab'
%   File extension for data files found in ``fullfile('MAG', sc)``, beginning with ``'.'``.
% SEQUENTIAL : bool
%   Whether to plot points by index or hours relative to a reference time.
% coeffPath : char, 1xD, default='modelCoeffs'
%   Directory containing model coefficients files.
% figDir : char, 1xE, default='figures'
%   Directory to use for output figures.
% figXtn : char, 1xF, default='pdf'
%   Extension to use for figures, which determines the file type.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('scName', 'var'); scName = "Swarm"; end
    if ~exist('xtn', 'var'); xtn = '.tab'; end
    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 1; end
    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    
    parentName = 'Earth';
    sc = char(scName);
    SetPlotDefaults();
    
    cspice_kclear;
    fullOrbFormatSpec = '%23s%20f%20f%20f%20f%20f%20f%[^\n\r]';
    disp(['Importing data files over ' parentName ' orbits.'])
    orbStr = [parentName ' data'];
    files = dir(fullfile('MAG', sc, ['*' xtn]));
    % Swarm low-res (LR) decimated data are 1-day files with 1 s cadence, and there are 3
    % spacecraft
    nTot = 86400*3;
    [BrSC, BthSC, BphiSC, r_km, theta, phi] = deal(zeros(1, nTot));
    for iFile=1:length(files)
        datFile = files(iFile).name;
        disp(['Loading ' sc ' MAG data from ' datFile '.'])
        fileID = fopen(fullfile('MAG', sc, datFile), 'r');
        magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', ...
            'EndOfLine', '\r\n');
        fclose(fileID);
    
        iStart = 1 + 86400*(iFile-1);
        iEnd = iStart + 86399;
        t_UTC(iStart:iEnd) = magData{1}';
        BrSC(iStart:iEnd) = magData{2}';
        BthSC(iStart:iEnd) = magData{3}';
        BphiSC(iStart:iEnd) = magData{4}';
        r_km(iStart:iEnd) = magData{5}';
        theta(iStart:iEnd) = magData{6}';
        phi(iStart:iEnd) = magData{7}';
    end
    
    xyz_km = zeros(3, nTot);
    xyz_km(1,:) = r_km .* sin(theta) .* cos(phi);
    xyz_km(2,:) = r_km .* sin(theta) .* sin(phi);
    xyz_km(3,:) = r_km .* cos(theta);
    
    spkParent = LoadSpice(parentName, sc);
    disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
    ets = cspice_str2et(t_UTC);
    t_h = ets / 3600;
    
    % Plot latitude to identify correlations with field model deviations
    if SEQUENTIAL
        xx = 1:nTot;
        xInfo = 'Measurement index';
    else
        xx = t_h - 175308.0192178;
        xInfo = 'Time relative to NY 2020 (h)';
    end
    lat_deg = 90 - rad2deg(theta);

    windowName = [char(scName) ' latitudes'];
    yy = lat_deg;
    yInfo = 'Latitude (degrees)';
    titleInfo = 'Swarm ABC latitudes on 2020-01-01';
    legendStrings = "latitude";
    fName = 'SwarmABClatitudes';
    trajFnum = 9001;
    PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, trajFnum);
    
    %% Plot and calculate products
    nOpts = 1; nMPopts = 0;
    opts = 1:nOpts;
    MPopts = -1:-1;
    for opt=opts
        for MPopt=MPopts
            PlotBandLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, scName, ...
                parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL, coeffPath, figDir, ...
                figXtn, LIVE_PLOTS);
        end
    end

end