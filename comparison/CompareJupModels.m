function CompareJupModels(LIVE_PLOTS, scName, moonName, orbNum, moonProx_RP, PlanetMinDist_RP, ...
    PlanetMaxDist_RP, dataDir, figDir, figXtn, SEQUENTIAL, coeffPath, FULLORBITS, FLYBYS, JUNOTOO)
% Compare magnetic field measurements from spacecraft near Jupiter against each implemented
% magnetic field model.
%
% Datasets:
%
%   * Galileo MAG: https://doi.org/10.17189/1519667, volume GO-J-MAG-3-RDR-HIGHRES-V1.0
%   * Juno MAG: https://doi.org/10.17189/1519711, volume JNO-J-3-FGM-CAL-V1.0
%
% Note
% ----
% This script requires having run PlanetMag.m for the Galileo era for all implemented planetary and
% magnetopause models using ``ALL_MODELS=1``.
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% scName : string, default="Galileo"
%   Spacecraft name for which magnetic field data will be compared against implemented models. A
%   directory must exist with this name in the ``MAG`` directory within the top-level P\lanetMag
%   directory. This directory will be searched for body-specific data files, and each of these
%   files will be loaded.
% moonName : char, 1xC, default='Ganymede'
%   Name of moon for which ``scName`` has flyby data. Flyby data is not used if ``FULLORBITS`` is
%   true, but a name must be passed for kernel loading.
% orbNum : int, default=-1
%   Orbit number to use for data comparison. Accepts the following special options:
%
%   - ``-1``: Use all orbits for which data is available.
%   - ``-2``: Use all orbits for which flyby data is available for ``moonName``.
%
% moonProx_RP : double, default=0.1
%   Minimum distance to major moons of the system within which measurements will be ignored in
%   model validation due to local perturbations, in units of the parent planet's radius.
% PlanetMinDist_RP : double, default=10
%   Minimum distance to parent planet within which measurements will be ignored, in units of the
%   parent planet's radius.
% PlanetMaxDist_RP : double, default=20
%   Maximum distance away from parent planet beyond which measurements will be ignored, in units of
%   the parent planet's radius.
% dataDir : char, 1xD, default='out'
%   Name of directory where excitation moments are printed to disk.
% figDir : char, 1xE, default='figures'
%   Directory to use for output figures.
% figXtn : char, 1xF, default='pdf'
%   Extension to use for figures, which determines the file type.
% SEQUENTIAL : bool, default=1
%   Whether to plot points by index or hours relative to a reference time (J2000).
% coeffPath : char, 1xG, default='modelCoeffs'
%   Directory containing model coefficients files.
% FULLORBITS : bool, default=1
%   Whether to evalate goodness of fit from full-orbit data in System III coordinates beyond some
%   threshold distance ``rMinMoon_RP``. Independent of ``FLYBYS``.
% FLYBYS : bool, default=0
%   Whether to evaluate goodness of fit from flyby data in moon coordinates. Measurements are only
%   used if they are beyond ``rMinMoon_RP`` away from the target moon. Independent of
%   ``FULLORBITS``.
% JUNOTOO : bool, default=1
%   Whether to include Juno flyby data if ``FLYBYS`` is also true. Currently, only Ganymede flyby
%   data is implemented.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('scName', 'var'); scName = "Galileo"; end
    if ~exist('moonName', 'var'); moonName = 'Ganymede'; end
    if ~exist('orbNum', 'var'); orbNum = -1; end
    if ~exist('moonProx_RP', 'var'); moonProx_RP = 0.1; end
    if ~exist('PlanetMinDist_RP', 'var'); PlanetMinDist_RP = 10; end
    if ~exist('PlanetMaxDist_RP', 'var'); PlanetMaxDist_RP = 20; end
    if ~exist('dataDir', 'var'); dataDir = 'out'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 1; end
    if ~exist('FULLORBITS', 'var'); FULLORBITS = 1; end
    if ~exist('FLYBYS', 'var'); FLYBYS = 0; end
    if ~exist('JUNOTOO', 'var'); JUNOTOO = 1; end

    cspice_kclear;
    parentName = 'Jupiter';
    
    LoadSpice(moonName, char(scName));
    if strcmp(scName, "Galileo")
        sc = 'Galileo Orbiter';
    else
        sc = char(scName);
    end
    
    switch moonName
        case 'Io'
            fbList = [0, 24, 27, 31, 32, 33];
            fbCode = 'IO';
        case 'Europa'
            fbList = [4, 11, 12, 14, 15, 19, 26];
            fbCode = 'EUR';
        case 'Ganymede'
            fbList = [1, 2, 7, 8, 28, 29];
            fbCode = 'GAN';
        case 'Callisto'
            fbList = [3, 9, 10, 30];
            fbCode = 'CALL';
    end
    
    if orbNum == -1
        orbStr = 'all orbits';
        fbStr = 'all flybys';
        orbList = linspace(0, 35, 36);
        orbList(orbList == 5) = []; % Solar conjunction on orbit 5 means no data at all
    elseif orbNum == -2
        orbStr = ['all orbits with ' moonName ' flybys'];
        fbStr = 'all flybys';
        orbList = fbList;
    else
        orbStr = ['orbit ' num2str(orbNum)];
        fbStr = ['flyby ' moonName(1) num2str(orbNum)];
        orbList = orbNum;
        fbList = orbNum;
    end
    
    if FULLORBITS && strcmp(scName, "Galileo")
        [t_UTC, BrSC, BthSC, BphiSC] = deal([]);
        fullOrbFormatSpec = '%23s%24s%10f%10f%10f%10f%7f%7f%7f%f%[^\n\r]';
        disp(['Importing PDS files over ' orbStr '.'])
        for i=1:length(orbList)
            if orbNum ~= -1
                disp(['Loading ' char(scName) ' MAG data for orbit ' num2str(orbList(i)) '.'])
            end
            datFile = fullfile(['MAG/' char(scName) '/Jupiter/ORB' sprintf('%02d', orbList(i)) ...
                '_SYS3.TAB']);
            fileID = fopen(datFile,'r');
            magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', ...
                'char', 'EndOfLine', '\r\n');
            fclose(fileID);
    
            t_UTC = [t_UTC, magData{1}'];
            BrSC = [BrSC, magData{3}'];
            BthSC = [BthSC, magData{4}'];
            BphiSC = [BphiSC, magData{5}'];
        end
    
        LoadSpice(moonName, char(scName));
        R_P = 71492;
    
        nTot = length(t_UTC);
        disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
        ets = cspice_str2et(t_UTC);
        disp(['Getting moon- and Jupiter-relative distances for all ' num2str(nTot) ' points.'])
        [rMinMoon_km, rJup_km] = GetMinMoonDist(sc, parentName, ets);
    
        % Delete measurement times near moons, far from Jupiter, and junk data
        BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
        finiteMax_nT = 2500;
        RPunit = ' R_J';
        disp(['Excluding all points satisfying at least one of the following:' newline ...
              'Distance to a major moon < ' num2str(moonProx_RP) RPunit newline ...
              'Planetocentric distance < ' num2str(PlanetMinDist_RP) RPunit newline ...
              'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
              'Suspect measurements, |B| > ' num2str(finiteMax_nT) 'nT.'])
        % Full limits
        exclude = find(rMinMoon_km/R_P < moonProx_RP | rJup_km/R_P < PlanetMinDist_RP ...
            | rJup_km/R_P > PlanetMaxDist_RP | BmagSC > finiteMax_nT);
        % No parent planet distance limits (full magnetosphere)
        %exclude = find(rMinMoon_km/R_P < moonProx_RP | BmagSC > finiteMax_nT);
        ets(exclude) = [];
        BrSC(exclude) = [];
        BthSC(exclude) = [];
        BphiSC(exclude) = [];
    
        npts = length(ets);
        t_h = ets / 3600;
        disp(['Getting ' sc ' positions for ' num2str(npts) ' pts over ' orbStr '.'])
        [r_km, theta, phi, xyz_km, spkParent] = GetPosSpice(sc, parentName, t_h);
    else
        spkParent = ['IAU_' upper(parentName)];
    end
    spkMoon = ['IAU_' upper(moonName)];
    
    %% Flyby data
    
    if FLYBYS
        [fbt_UTC, fbBrSC, fbBthSC, fbBphiSC] = deal([]);
        flybyFormatSpec = '%23s%10f%10f%10f%10f%7f%7f%7f%f%[^\n\r]';
        for i=1:length(fbList)
            if orbNum ~= -1
                disp(['Loading ' char(scName) ' MAG data for flyby ' moonName(1) ...
                    num2str(orbList(i)) '.'])
            end
            datFile = fullfile(['MAG/' char(scName) '/' moonName '/ORB' sprintf('%02d_', ...
                fbList(i)) fbCode '_SYS3.TAB']);
            fileID = fopen(datFile,'r');
            magData = textscan(fileID, flybyFormatSpec, Inf, 'Delimiter', '', 'TextType', ...
                'char', 'EndOfLine', '\r\n');
            fclose(fileID);
    
            fbt_UTC = [fbt_UTC, magData{1}'];
            fbBrSC = [fbBrSC, magData{2}'];
            fbBthSC = [fbBthSC, magData{3}'];
            fbBphiSC = [fbBphiSC, magData{4}'];
        end
    
        fbets = cspice_str2et(fbt_UTC);
        fbt_h = fbets / 3600;
        [fbr_km, fbtheta, fbphi, fbxyz_km, ~] = GetPosSpice(sc, parentName, fbt_h);
        
        [BxSCS3, BySCS3, BzSCS3] = Bsph2Bxyz(fbBrSC, fbBthSC, fbBphiSC, fbtheta, fbphi);
        [BxSC, BySC, BzSC] = RotateVecSpice(BxSCS3, BySCS3, BzSCS3, fbets, spkParent, spkMoon);
        r_RM = GetTargetMoonDist(sc, moonName, parentName, fbets);
    end
    
    % Add Juno flyby data
    if JUNOTOO && FLYBYS
        LoadSpice(moonName, 'Juno');
        fbjNum = 34;
        datFile = fullfile(['MAG/Juno/' moonName '/ORB' num2str(fbjNum) '_FGM_IAU.txt']);
        disp(['Loading FGM data for Juno ' moonName ' flyby on orbit ' num2str(fbjNum) '.'])
        junoFlybyFormatSpec = '%23s%f%f%f%f%f%f%f%[^\n\r]';
        fileID = fopen(datFile,'r');
        magData = textscan(fileID, junoFlybyFormatSpec, Inf, 'Delimiter', '', 'TextType', ...
            'char', 'EndOfLine', '\r\n');
        
        jt_UTC = magData{1}';
        jBx = magData{2}';
        jBy = magData{3}';
        jBz = magData{4}';
        jets = cspice_str2et(jt_UTC);
        jt_h = jets / 3600;
        
        [jr_km, jtheta, jphi, jxyz_km, ~] = GetPosSpice('Juno', parentName, jt_h);
        jr_RM = GetTargetMoonDist('Juno', moonName, parentName, jets);
        
        % Append to lists of coords to eval in Jupiter field model
        fbets = [fbets, jets];
        % NOT t_h, we will use the separate arrays to organize outputs.
        fbr_km = [fbr_km, jr_km];
        fbtheta = [fbtheta, jtheta];
        fbphi = [fbphi, jphi];
        fbxyz_km = [fbxyz_km, jxyz_km];
        r_RM = [r_RM, jr_RM];
        BxSC = [BxSC, jBx];
        BySC = [BySC, jBy];
        BzSC = [BzSC, jBz];
        
        scName = [scName, "Juno"];
    else
        jt_h = [];
    end
    
    %% Plot and calculate products
    nOpts = 7; nMPopts = 5;
    opts = 1:nOpts;
    MPopts = 1:(nMPopts + 1); % Add 1 to force noMP model in addition
    for opt=opts
        for MPopt=MPopts
            if FULLORBITS
                PlotBandLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, scName, ...
                    parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL, coeffPath, figDir, ...
                    figXtn, LIVE_PLOTS);
            else
                PlotBandLsqMoon(fbets, fbt_h, fbr_km, fbtheta, fbphi, fbxyz_km, r_RM, BxSC, ...
                    BySC, BzSC, scName, parentName, spkParent, moonName, char(scName), fbStr, ...
                    opt, MPopt, SEQUENTIAL, dataDir, figDir, figXtn, LIVE_PLOTS, jt_h);
            end
        end
    end

end