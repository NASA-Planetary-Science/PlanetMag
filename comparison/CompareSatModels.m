function CompareSatModels(LIVE_PLOTS, scNames, SEQUENTIAL, coeffPath, figDir, figXtn)
% Compare magnetic field measurements from spacecraft near Saturn against each implemented magnetic
% field model.
%
% Datasets:
%
%   * Voyager 1 MAG: https://doi.org/10.17189/1522661, volume VG1-S-MAG-4-SUMM-L1COORDS-48SEC-V1.0
%   * Voyager 2 MAG: https://doi.org/10.17189/1520002, volume VG2-S-MAG-4-SUMM-L1COORDS-48SEC-V1.1
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% scNames : string, 1xS, default=["Voyager 1", "Voyager 2"]
%   Spacecraft names for which magnetic field data will be compared against implemented models. A
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

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('scNames', 'var'); scNames = ["Voyager 1", "Voyager 2"]; end
    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 0; end
    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    
    cspice_kclear;
    parentName = 'Saturn';
    SetPlotDefaults();
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt
    
    fullOrbFormatSpec = ['%23s%12s%1d%9f%9f%9f%9f%9f%7f%7f%8f%8f%8f%10s%[^\n\r]', ...
        '%24s%6d%3d%2d%2d%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%2d%10s%[^\n\r]'];
    disp(['Importing PDS files over ' parentName ' flybys.'])
    orbStr = [parentName ' flyby'];
    Rp_km = 60268;
    for i = 1:2
        LoadSpice(parentName, char(scNames(i)));
    end
    etStime = ['1980-11-12T23:45:42.725'; '1981-08-26T03:24:04.760'];
    etS = cspice_str2et(etStime);
    tCA_h = [-167724.2239, -160856.5842];
    
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

    the = linspace(0,pi,50); ph = linspace(0,2*pi,100);
    [the2D, ph2D] = meshgrid(the,ph);
    xp = sin(the2D) .* cos(ph2D); yp = sin(the2D) .* sin(ph2D); zp = cos(the2D);
    surf(xp,yp,zp, 'FaceColor', '#EDB120');
    axlim = 30;
    pbaspect([1 1 1]);
    xlim([-axlim, axlim]); ylim([-axlim, axlim]); zlim([-axlim, axlim]);
    grid on;
    title([parentName ' flyby trajectories']);
    xlabel('x KSO (R_S, toward Sun)');
    ylabel('y KSO (R_S), \approx orbital v');
    zlabel('z KSO (R_S)');
    
    rotMat = cspice_pxform('IAU_SATURN', 'KSO', etS(1));
    vecMat = zeros(3,1);
    vecMat(3,1) = 2.5;
    poleKSO = squeeze(pagemtimes(rotMat, vecMat));
    plot3([0,poleKSO(1)], [0,poleKSO(2)], [0,poleKSO(3)], 'Color', 'r', 'LineWidth', 2)
    scatter3(poleKSO(1), poleKSO(2), poleKSO(3), 15, 'r')
    
    xi = 0.72;
    Rss = 19;
    thDSZ = linspace(0,pi,101);
    thDSZ = thDSZ(1:end-1);
    phMP = linspace(0,2*pi,100);
    [th2D, ph2D] = meshgrid(thDSZ, phMP);
    rMP = Rss .* (2 ./ (1 + cos(th2D))).^xi;
    xMP = rMP .* cos(thDSZ) .* ones(size(ph2D));
    yMP = rMP .* sin(thDSZ) .* sin(ph2D);
    zMP = rMP .* sin(thDSZ) .* cos(ph2D);
    MPsurf = [xMP; yMP; zMP];
    surf(xMP, yMP, zMP, 'FaceColor', 'b')
    
    for i=2:-1:1
        sc = char(scNames(i));
        datFile = fullfile(['MAG/' sc '/vg' num2str(i) '_' parentName '_48s_sph.tab']);
        disp(['Loading ' sc ' MAG data from ' datFile(i) '.'])
        fileID = fopen(datFile,'r');
        magData = textscan(fileID, fullOrbFormatSpec, Inf, 'Delimiter', ',', 'TextType', ...
            'char', 'EndOfLine', '\r\n');
        fclose(fileID);
    
        t_UTC{i} = magData{1}';
        if i == 1
            BrSC{i} = magData{4}';
            BthSC{i} = magData{5}';
            BphiSC{i} = magData{6}';
        else
            BrSC{i} = magData{6}';
            BthSC{i} = magData{7}';
            BphiSC{i} = magData{8}';
        end
    
        nTot = length(t_UTC{i});
        disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
        ets{i} = cspice_str2et(t_UTC{i});
        disp(['Getting moon and planet distances for all ' num2str(nTot) ' points.'])
        [rMinMoon_km, rP_km{i}] = GetMinMoonDist(sc, parentName, ets{i});
    
        % Delete measurement times far from Saturn and junk data
        BmagSC{i} = sqrt(BrSC{i}.^2 + BthSC{i}.^2 + BphiSC{i}.^2);
        moonProx_RP = 0.1;
        PlanetMaxDist_RP = 60;
        finiteMax_nT = 8e3;
        RPunit = [' R_' parentName(1)];
        disp(['Excluding all points satisfying at least one of the following:' newline ...
            'Distance to a major moon < ' num2str(moonProx_RP) RPunit newline ...
            'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
            'Suspect measurements, |B| > ' num2str(finiteMax_nT) 'nT.'])
        % Full limits
        exclude = find(rMinMoon_km/Rp_km < moonProx_RP | rP_km{i}/Rp_km > PlanetMaxDist_RP ...
            | BmagSC{i} > finiteMax_nT);
        ets{i}(exclude) = [];
        BrSC{i}(exclude) = [];
        BthSC{i}(exclude) = [];
        BphiSC{i}(exclude) = [];
    
        npts = length(ets{i});
        t_h{i} = ets{i} / 3600;
        disp(['Getting ' sc ' positions for ' num2str(npts) ' pts.'])
        [r_km, theta, phi, xyz_km, spkParent] = GetPosSpice(sc, parentName, t_h{i});
    
        %% Plot and calculate products
        nOpts = 2; nMPopts = 2;
        opts = 0:0;
        MPopts = 2:(nMPopts + 1); % Add 1 to force noMP model in addition
        for opt=opts
            for MPopt=MPopts(2:end)
                PlotBandLsq(ets{i}, t_h{i}, r_km, theta, phi, xyz_km, BrSC{i}, BthSC{i}, ...
                    BphiSC{i}, scNames(i), parentName, spkParent, orbStr, opt, MPopt, ...
                    SEQUENTIAL, coeffPath, figDir, figXtn, LIVE_PLOTS);
            end
        end
    
        %% Plot trajectory
        % t_h = linspace(cspice_str2et('1986-01-23T00:00:00.000'), ...
        %     cspice_str2et('1986-01-27T00:00:00.000'), 5000) / 3600;
        [~, ~, ~, xyzKSO_km, ~] = GetPosSpice(sc, parentName, t_h{i}, 'KSO');
        xyz_Rp = xyzKSO_km / Rp_km;
        x = xyz_Rp(1,:); y = xyz_Rp(2,:); z = xyz_Rp(3,:);
        % Re-focus existing trajectory plot
        if LIVE_PLOTS
            fig = figure(trajFnum);
            set(gcf, 'Visible', 'on', 'Name', windowName);
        else
            fig = figure('Visible', 'off', 'Name', windowName);
        end
        hold on;
        scTraj{i} = plot3(x,y,z, 'LineWidth', 1.5, 'DisplayName', sc);
        
    end
    legend([scTraj{1} scTraj{2}], scNames)

    if ~LIVE_PLOTS
        outFig = fullfile(figDir, [parentName 'Trajectories.' figXtn]);
        fig.Units = fig.PaperUnits;
        fig.PaperSize = fig.Position(3:4);
        saveas(fig, outFig)
        disp(['Figure saved to ' outFig '.'])
    end
    close(fig)
end