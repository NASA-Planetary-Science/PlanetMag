function PlotBandLsqMoon(ets, t_h, r_km, theta, phi, xyz_km, r_RM, BxSC, BySC, BzSC, scName, ...
    parentName, S3coords, moonName, era, fbStr, opt, MPopt, SEQUENTIAL, dataDir, figDir, ...
    figXtn, LIVE_PLOTS, coeffPath, jt_h)
% Plots and calculates comparisons between modeled and measured magnetic fields in the vicinity of
% a target moon.
%
% Generates time series data of a specified combination of magnetic field models implemented in
% P\lanetMag and compares against spacecraft measurements of the planetary magnetic field.
% Comparisons are plotted and least-squares differences are calculated and printed to the terminal.
% Intended as a final step in model validation; unlike PlotBandLsq, this function uses flyby data
% in the vicinity of the target moon and the recorded excitation moments to model the magnetic
% field applied by the parent planet.
%
% Note
% ----
% Because this function uses excitation moments saved from PlanetMag, that function must be run for
% each of the input settings (``moonName``, ``opt``, ``MPopt``, and ``era``) before this one.
%
% Parameters
% ----------
% ets : double, 1xN
%   Ephemeris times of measurements to compare in TDB seconds relative to J2000.
% t_h : double, 1xN
%   Ephemeris times of measurements to compare in TDB hours relative to J2000 (``ets/3600``).
% r_km : double, 1xN
%   Radial distance of measurement locations from planet center of mass in km.
% theta : double, 1xN
%   Colatitude of measurement locations from planet center of mass in radians.
% phi : double, 1xN
%   East longitude of measurement locations from planet center of mass in radians.
% xyz_km : double, 3xN
%   Cartesian coordinates of measurement locations in ``S3coords`` frame in km.
% r_RM : double, 1xN
%   Radial distance of measurement locations from moon center of mass in moon equatorial radii.
% BxSC : double, 1xN
%   x component (:math:`B_x`) of magnetic field measurements in ``S3coords`` frame in nT.
% BySC : double, 1xN
%   y component (:math:`B_y`) of magnetic field measurements in ``S3coords`` frame in nT.
% BzSC : double, 1xN
%   z component (:math:`B_z`) of magnetic field measurements in ``S3coords`` frame in nT.
% scName : string, 1xS'
%   Name(s) of spacecraft for measurement comparisons. Accepts a lone string or a list of strings.
% parentName : char, 1xC
%   Name of parent planet the evaluated model(s) and spacecraft measurements apply to.
% S3coords : char, 1xD
%   Standard coordinate frame (for the parent planet) used for evaluation of magnetic fields.
% moonName : char, 1xE
%   Name of target moon the evaluated model(s) and spacecraft measurements apply to.
% era : char, 1xF
%   Name of spacecraft mission time frame to use for excitation moments. Must be 
% fbStr : char, 1xG
%   Description of flyby(s) covered to place in legend labels.
% opt : int
%   Index of planetary magnetic field model. See GetModelOpts for more details and available
%   options.
% MPopt : int
%   Index of magnetopause current magnetic field model. See GetModelOpts for more details and
%   available options.
% SEQUENTIAL : bool, default=0
%   Whether to plot points by index or hours relative to a reference time (closest approach).
% dataDir : char, 1xH, default='out'
%   Name of directory where excitation moments are printed to disk.
% figDir : char, 1xI, default='figures'
%   Directory to use for output figures.
% figXtn : char, 1xF, default='pdf'
%   Extension to use for figures, which determines the file type.
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% coeffPath : char, 1xJ, default='modelCoeffs'
%   Directory containing model coefficients files.
% jt_h : double, 1xM, default=[]
%   Ephemeris times of Juno measurements in TDB hours relative to J2000 to use in data comparisons.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('dataDir', 'var'); dataDir = 'out'; end
    if ~exist('coeffpath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('jt_h', 'var'); jt_h = []; end
    if isempty(jt_h); JUNOTOO=0; else; JUNOTOO=1; end
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt
    
    [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt, MPopt);
    magPhase = 0;
    npts = length(ets);
    
    disp(['Evaluating ' magModelDescrip ' for flybys.'])
    if strcmp(magModelDescrip, 'KS2005')
        [Bvec, Mdip_nT, Odip_km] = MagFldJupiterKS2005(r_km, theta, phi, ets, 1);
    else
        [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, r_km, theta, phi, MagModel, ...
            CsheetModel, magPhase, 1);
    end
    if ~(strcmp(MPmodel, 'None') || contains(magModelDescrip, 'KS2005'))

        nSW_pcc = 0.14 * ones(1,npts);
        vSW_kms = 400  * ones(1,npts);
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, [t_h, jt_h]*3600, xyz_km, Mdip_nT, ...
            Odip_km, parentName, S3coords, MPmodel, coeffPath, 1);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

    end
    
    spkMoon = ['IAU_' upper(moonName)];
    [BxS3, ByS3, BzS3] = Bsph2Bxyz(Bvec(1,:), Bvec(2,:), Bvec(3,:), theta, phi);
    [Bx, By, Bz] = RotateVecSpice(BxS3, ByS3, BzS3, ets, S3coords, spkMoon);
    
    % Get excitation moments from longer time series with each model
    excMomentsFile = fullfile(dataDir, ['Be1xyz_' moonName '_' era '_' fEnd '.txt']);
    reloadData = dlmread(excMomentsFile, ',', 1, 0);

    Texc_h = reloadData(:,1);
    B0x = mean(reloadData(:,2));
    B0y = mean(reloadData(:,3));
    B0z = mean(reloadData(:,4));
    B1x_Re = reloadData(:,5);
    B1x_Im = reloadData(:,6);
    B1y_Re = reloadData(:,7);
    B1y_Im = reloadData(:,8);
    B1z_Re = reloadData(:,9);
    B1z_Im = reloadData(:,10);
    npeaks = length(Texc_h);
    
    omega_ph = 2*pi./Texc_h;
    B1x = B1x_Re + 1i*B1x_Im;
    B1y = B1y_Re + 1i*B1y_Im;
    B1z = B1z_Re + 1i*B1z_Im;

    npts_noJuno = length(t_h);
    [BxExc, ByExc, BzExc] = deal(ones(1,npts_noJuno));
    BxExc = B0x * BxExc;
    ByExc = B0y * ByExc;
    BzExc = B0z * BzExc;
    
    for i=1:npeaks
        BxExc = BxExc + real(B1x(i) * exp(-1i * omega_ph(i) * t_h));
        ByExc = ByExc + real(B1y(i) * exp(-1i * omega_ph(i) * t_h));
        BzExc = BzExc + real(B1z(i) * exp(-1i * omega_ph(i) * t_h));
    end
    
    %% Add Juno data if applicable    
    if JUNOTOO
        % Get excitation moments from longer time series with each model
        excMomentsFile = fullfile(dataDir, ['Be1xyz_' moonName '_Juno_' fEnd '.txt']);
        reloadData = dlmread(excMomentsFile, ',', 1, 0);

        Texc_h = reloadData(:,1);
        B0x = mean(reloadData(:,2));
        B0y = mean(reloadData(:,3));
        B0z = mean(reloadData(:,4));
        B1x_Re = reloadData(:,5);
        B1x_Im = reloadData(:,6);
        B1y_Re = reloadData(:,7);
        B1y_Im = reloadData(:,8);
        B1z_Re = reloadData(:,9);
        B1z_Im = reloadData(:,10);
        npeaks = length(Texc_h);

        omega_ph = 2*pi./Texc_h;
        B1x = B1x_Re + 1i*B1x_Im;
        B1y = B1y_Re + 1i*B1y_Im;
        B1z = B1z_Re + 1i*B1z_Im;

        nptsJuno = length(jt_h);
        [jBxExc, jByExc, jBzExc] = deal(ones(1,nptsJuno));
        jBxExc = B0x * jBxExc;
        jByExc = B0y * jByExc;
        jBzExc = B0z * jBzExc;

        for i=1:npeaks
            jBxExc = jBxExc + real(B1x(i) * exp(-1i * omega_ph(i) * jt_h));
            jByExc = jByExc + real(B1y(i) * exp(-1i * omega_ph(i) * jt_h));
            jBzExc = jBzExc + real(B1z(i) * exp(-1i * omega_ph(i) * jt_h));
        end
        
        BxExc = [BxExc, jBxExc];
        ByExc = [ByExc, jByExc];
        BzExc = [BzExc, jBzExc];
        
        % Now combine the t_h arrays for plotting
        t_h = [t_h, jt_h];
    end
    
    %%
    
    excStr = "Excitation moments";
    commonTitle = [parentName ' model comparison vs ' char(scName) ' data'];
    npts = length(t_h);
    if SEQUENTIAL
        xx = 1:npts;
        xInfo = 'Measurement index';
    else
        xx = t_h;
        xInfo = 'Time past J2000 (h)';
    end
    if JUNOTOO; xJuno = [xx(npts_noJuno+1), xx(npts)]; end
    
    figNumBase = 3000 + 100*opt + 10*MPopt;
    windowName = ['Bx, ' fbStr ', ' magModelDescrip];
    yy = [Bx; BxSC; BxExc];
    yInfo = 'Vector component (nT)';
    titleInfo = commonTitle;
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    legendStrings = [string(magModelDescrip), strcat(scName, ' MAG'), excStr];
    if JUNOTOO
        vlinexCols = {xJuno, '#B000FF'};
    else
        vlinexCols = {};
    end
    fName = [moonName 'BxFlybyComparison' magModelDescrip];
    fig = PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 1, 'linear', 'linear', 'auto', 'auto', {}, 0, ...
        0, vlinexCols);
    close(fig);

    windowName = ['By, ' fbStr ', ' magModelDescrip];
    yy = [By; BySC; ByExc];
    fName = [moonName 'ByFlybyComparison' magModelDescrip];
    fig = PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 2, 'linear', 'linear', 'auto', 'auto', {}, 0, ...
        0, vlinexCols);
    close(fig);

    windowName = ['Bz, ' fbStr ', ' magModelDescrip];
    yy = [Bz; BzSC; BzExc];
    fName = [moonName 'BzFlybyComparison' magModelDescrip];
    fig = PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 3, 'linear', 'linear', 'auto', 'auto', {}, 0, ...
        0, vlinexCols);
    close(fig);

    BxD = Bx - BxSC;
    ByD = By - BySC;
    BzD = Bz - BzSC;
    BxDexc = BxExc - BxSC;
    ByDexc = ByExc - BySC;
    BzDexc = BzExc - BzSC;

    windowName = ['IAU vector comp diffs, ' fbStr ', ' magModelDescrip ' - MAG'];
    yy = [BxD; ByD; BzD; BxDexc; ByDexc; BzDexc];
    yInfo = 'Component difference (nT)';
    legendStrings = [string([mathTxt '\Delta B_x']), string([mathTxt '\Delta B_y']), ...
        string([mathTxt '\Delta B_z']), string([mathTxt '\Delta B^e_x']), string([mathTxt ...
        '\Delta B^e_y']), string([mathTxt '\Delta B^e_z'])];
    fName = [moonName 'DeltaBflyby' magModelDescrip];
    fig = PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 4, 'linear', 'linear', 'auto', 'auto', {}, 0, ...
        0, vlinexCols);
    close(fig);

    RMmin = 2;
    iFar = find(r_RM > RMmin);
    nptsFar = length(iFar);
    BxLsq = BxD(iFar).^2;
    ByLsq = ByD(iFar).^2;
    BzLsq = BzD(iFar).^2;
    chi2 = sum([BxLsq, ByLsq, BzLsq], 'all') / 3 / nptsFar;
    disp(['Across flybys, beyond ' num2str(RMmin) ...
        ' R_' moonName(1) ', for this model chi^2 = ' sprintf('%.2f', chi2) '.'])
    BxExcLsq = BxDexc(iFar).^2;
    ByExcLsq = ByDexc(iFar).^2;
    BzExcLsq = BzDexc(iFar).^2;
    chi2exc = sum([BxExcLsq, ByExcLsq, BzExcLsq], 'all') / 3 / nptsFar;
    disp(['Excitation moments only, beyond ' num2str(RMmin) ...
        ' R_' moonName(1) ': chi^2 = ' sprintf('%.2f', chi2exc) '.'])
end
