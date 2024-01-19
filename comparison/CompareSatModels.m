function CompareSatModels(LIVE_PLOTS, LOAD_PDS_ASCII, yearRange, RELATIVE_t, RELATIVE_r, ...
    scName, SEQUENTIAL, coeffPath, figDir, figXtn, magDataMat)
% Compare magnetic field measurements from spacecraft near Saturn against each implemented magnetic
% field model.
%
% Datasets:
%
%   * Cassini MAG: 	https://doi.org/10.17189/5rhj-sm88, volume CO-E/SW/J/S-MAG-4-SUMM-1MINAVG-V2.1
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% LOAD_PDS_ASCII : bool, default=0
%   Whether to load in .TAB files downloaded from PDS (true) or use pre-converted and compressed
%   .mat summery of just the important bits (false).
% yearRange : int, 1xG, default=4:17
%   Years over which to compare Cassini data. Default includes all measurements.
% RELATIVE_t : bool, default=0
%   Whether to plot points relative to the start of the input time series. Only has an effect when
%   ``SEQUENTIAL = 0``.
% RELATIVE_r : bool, default=0
%   Whether to plot points relative to distance from the parent planet. Overrides SEQUENTIAL and
%   RELATIVE_t.
% scName : string, 1xS, default=["Cassini"]
%   Spacecraft name for which magnetic field data will be compared against implemented models. A
%   directory must exist with each name in the ``MAG`` directory within the top-level P\lanetMag
%   directory. These directories will be searched for data files with the ``.tab`` extension, and
%   each of these files will be loaded.
% SEQUENTIAL : bool, default=0
%   Whether to plot points by index or hours relative to a reference time (typically closest
%   approach).
% coeffPath : char, 1xC, default='modelCoeffs'
%   Directory containing model coefficients files.
% figDir : char, 1xD, default='figures'
%   Directory to use for output figures.
% figXtn : char, 1xE, default='pdf'
%   Extension to use for figures, which determines the file type.
% magDataMat : char, 1xF, default='MAG/scName/ALL_FGM_KRTP_1M'
%   File path to use for save/reload of PDS data in .mat file.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('LOAD_PDS_ASCII', 'var'); LOAD_PDS_ASCII = 0; end
    if ~exist('yearRange', 'var'); yearRange = 4:17; end
    if ~exist('scName', 'var'); scName = "Cassini"; end
    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 0; end
    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    sc = char(scName);
    if ~exist('magDataMat', 'var'); magDataMat = fullfile('MAG', sc, 'ALL_FGM_KRTP_1M'); end
    if ~exist('RELATIVE_t', 'var'); RELATIVE_t = 0; end
    
    cspice_kclear;
    parentName = 'Saturn';
    SetPlotDefaults();
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt
    
    Rp_km = 60268;
    LoadSpice(parentName, sc);
    orbStr = [parentName ' orbit'];

    if LOAD_PDS_ASCII || ~exist([magDataMat '.mat'], 'file')

        magFormatSpec = '%19s%11f%11f%11f%11f%8f%7f%7f%5f%3d%[^\n\r]';
        disp(['Importing PDS measurements over ' parentName ' orbits.'])
        [t_UTC, BrSC, BthSC, BphiSC] = deal([]);
        for i=yearRange
            datFile = fullfile('MAG', sc, [num2str(2000+i) '_FGM_KRTP_1M.TAB']);
            disp(['Loading ' sc ' MAG data from ' datFile '.'])
            fileID = fopen(datFile,'r');
            magData = textscan(fileID, magFormatSpec, Inf, 'Delimiter', ',', 'TextType', ...
                'char', 'EndOfLine', '\r\n');
            fclose(fileID);
        
            t_UTC =  [t_UTC,  magData{1}'];
            BrSC =   [BrSC,   magData{2}'];
            BthSC =  [BthSC,  magData{3}'];
            BphiSC = [BphiSC, magData{4}'];
    
        end
        
        nTot = length(t_UTC);
        disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
        ets = cspice_str2et(t_UTC);
        save(magDataMat, 'ets', 'BrSC', 'BthSC', 'BphiSC');
        disp(['Saved PDS measurement times and field vectors to file: ' magDataMat])

    else
        disp(['Reloading PDS-archived measurement ephemeris times and field vectors from file:' ...
            ' ' magDataMat '.mat'])
        load(magDataMat, 'ets', 'BrSC', 'BthSC', 'BphiSC');
        nTot = length(ets);
    end

    disp(['Getting moon and planet distances for all ' num2str(nTot) ' points.'])
    [rMinMoon_km, rP_km] = GetMinMoonDist(sc, parentName, ets);

    % Delete measurement times far from Saturn, near moons, and junk data
    % Note: The B2010 model is derived only from measurements inside the Enceladus L shell,
    % i.e. r < 3.95 R_S.
    % Note: The Cassini 11 model intrinsic field is derived only from measurements in the FGM
    % range with |B| > 10e3 nT, and the current sheet model only from those in ranges with
    % |B| > 400 nT.
    BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
    moonProx_RP = 0.1;
    PlanetMaxDist_RP = 60;
    finiteMaxMag_nT = 20e3;
    RPunit = [' R_' parentName(1)];
    disp(['Excluding all points satisfying at least one of the following:' newline ...
        'Distance to a major moon < ' num2str(moonProx_RP) RPunit newline ...
        'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
        'Suspect measurements: |B| > ' num2str(finiteMaxMag_nT) 'nT,' newline ...
        '|B| > 400 when R_S > 4.50' 'nT.'])
    % Full limits
    rP_Rp = rP_km/Rp_km;
    exclude = find(rMinMoon_km/Rp_km < moonProx_RP | rP_Rp > PlanetMaxDist_RP ...
        | BmagSC > finiteMaxMag_nT | (BmagSC > 400 & rP_Rp > 4.50));
    ets(exclude) = [];
    BrSC(exclude) = [];
    BthSC(exclude) = [];
    BphiSC(exclude) = [];

    npts = length(ets);
    t_h = ets / 3600;
    disp(['Getting ' sc ' positions for ' num2str(npts) ' pts.'])
    [r_km, theta, phi, xyz_km, spkParent] = GetPosSpice(sc, parentName, t_h);

    %% Plot and calculate products
    nOpts = 2;
    opts = 1:nOpts;
    MPopts = 0:0; % Only noMP model
    for opt=opts
        for MPopt=MPopts
            PlotBandLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, scName, ...
                parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL, coeffPath, figDir, ...
                figXtn, LIVE_PLOTS, [], RELATIVE_t, RELATIVE_r);
        end
    end

    %% Plot trajectory
    disp(['Getting ' sc ' trajectories in KSO frame for ' num2str(npts) ' pts.'])
    [~, ~, ~, xyzKSO_km, ~] = GetPosSpice(sc, parentName, t_h, 'KSO');
    xyz_Rp = xyzKSO_km / Rp_km;
    x = xyz_Rp(1,:); y = xyz_Rp(2,:); z = xyz_Rp(3,:);
    
    %% Plot trajectories
    windowName = 'Spacecraft Saturn trajectories';
    trajFnum = 9001;
    if LIVE_PLOTS
        fig = figure(trajFnum);
        set(gcf, 'Visible', 'on', 'Name', windowName);
    else
        fig = figure('Visible', 'off', 'Name', windowName);
    end
    clf(); hold on;
    [interpreter, font] = SetPlotDefaults();
    ApplyPlotDefaults(fig, interpreter, font);
    
    % Plot planet, pole, and rings for illustrative purposes
    pbaspect([1 1 1]);
    axlim_Rp = 10;
    xlim([-axlim_Rp, axlim_Rp]); ylim([-axlim_Rp, axlim_Rp]); zlim([-axlim_Rp, axlim_Rp]);
    grid on;
    title([parentName ' orbit trajectories']);
    xlabel('x KSO (R_S, toward Sun)');
    ylabel('y KSO (R_S), \approx orbital v');
    zlabel('z KSO (R_S)');

    % Pole and IAU -> KSO rotation matrix
    rotMat = cspice_pxform('IAU_SATURN', 'KSO', ets(1));
    vecMat = zeros(3,1);
    vecMat(3,1) = 2.5;
    poleKSO = squeeze(pagemtimes(rotMat, vecMat));
    plot3([0,poleKSO(1)], [0,poleKSO(2)], [0,poleKSO(3)], 'Color', 'r', 'LineWidth', 2)
    scatter3(poleKSO(1), poleKSO(2), poleKSO(3), 15, 'r')

    % Plot planet surface for showing trajectories
    nth = 25; nph = 50;
    the = linspace(0,pi,nth); ph = linspace(0,2*pi,nph);
    [the2D, ph2D] = meshgrid(the,ph);
    xp = sin(the2D) .* cos(ph2D); yp = sin(the2D) .* sin(ph2D); zp = cos(the2D);
    xrlin = reshape(xp, 1, []); yrlin = reshape(yp, 1, []); zrlin = reshape(zp, 1, []);
    xyzSurf = [xrlin; yrlin; zrlin];
    xyzSurfRot = squeeze(pagemtimes(rotMat, xyzSurf));
    xp = reshape(xyzSurfRot(1,:), [nph, nth]);
    yp = reshape(xyzSurfRot(2,:), [nph, nth]);
    zp = reshape(xyzSurfRot(3,:), [nph, nth]);
    surf(xp,yp,zp, 'FaceColor', '#EDB120');

    % Plot rings for fun
    DringMin = 66970/Rp_km; % D ring min radius (all values here are from cpck31Oct2017.tpc)
    CringMin = 74510/Rp_km; % Bright, more visible C ring min radius
    AringMax = 136780/Rp_km; % Outer edge of largest visually contiguous portion
    FringMax = 140270/Rp_km; % Outer edge of outermost bright visible ring
    nr = 20; nph = 50;
    rRings = linspace(CringMin, AringMax, nr); ph = linspace(0,2*pi,nph);
    [r2D, ph2D] = meshgrid(rRings, ph);
    xr = r2D .* cos(ph2D); yr = r2D .* sin(ph2D); zr = zeros(size(ph2D));
    xrlin = reshape(xr, 1, []); yrlin = reshape(yr, 1, []); zrlin = reshape(zr, 1, []);
    xyzRings = [xrlin; yrlin; zrlin];
    xyzRKSO = squeeze(pagemtimes(rotMat, xyzRings));
    xRKSO = reshape(xyzRKSO(1,:), [nph, nr]);
    yRKSO = reshape(xyzRKSO(2,:), [nph, nr]);
    zRKSO = reshape(xyzRKSO(3,:), [nph, nr]);
    surf(xRKSO,yRKSO,zRKSO, 'FaceColor', '#6E6656');

    % Now plot trajectory
    scTraj = plot3(x,y,z, 'LineWidth', 1.5, 'DisplayName', sc);
    legend(scTraj, scName)

    outFig = fullfile(figDir, [parentName 'Trajectories.' figXtn]);
    fig.Units = fig.PaperUnits;
    fig.PaperSize = fig.Position(3:4);
    saveas(fig, outFig)
    disp(['Figure saved to ' outFig '.'])
    if ~LIVE_PLOTS; close(fig); end
end