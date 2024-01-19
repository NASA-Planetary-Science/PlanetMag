function PlotBandLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, scName, ...
    parentName, S3coords, orbStr, opt, MPopt, SEQUENTIAL, coeffPath, figDir, figXtn, ...
    LIVE_PLOTS, jt_h, RELATIVE_t, RELATIVE_r)
% Plots and calculates comparisons between modeled and measured magnetic fields.
%
% Generates time series data of a specified combination of magnetic field models implemented in
% P\lanetMag and compares against spacecraft measurements of the planetary magnetic field.
% Comparisons are plotted and least-squares differences are calculated and printed to the terminal.
% Intended as a final step in model validation.
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
% BrSC : double, 1xN
%   Radial component (:math:`B_r`) of magnetic field measurements in ``S3coords`` frame in nT.
% BthSC : double, 1xN
%   Colatitudinal component (:math:`B_\theta`) of magnetic field measurements in ``S3coords`` frame
%   in nT.
% BphiSC : double, 1xN
%   Azimuthal component (:math:`B_\phi`) of magnetic field measurements in ``S3coords`` frame in
%   nT.
% scName : string, 1xS'
%   Name(s) of spacecraft for measurement comparisons. Accepts a lone string or a list of strings.
% parentName : char, 1xC
%   Name of parent planet the evaluated model(s) and spacecraft measurements apply to.
% S3coords : char, 1xD
%   Standard coordinate frame used for evaluation of magnetic fields.
% orbStr : char, 1xE
%   Description of orbit(s) covered to place in legend labels.
% opt : int
%   Index of planetary magnetic field model. See GetModelOpts for more details and available
%   options.
% MPopt : int
%   Index of magnetopause current magnetic field model. See GetModelOpts for more details and
%   available options.
% SEQUENTIAL : bool, default=0
%   Whether to plot points by index or hours relative to a reference time.
% coeffPath : char, 1xE, default='modelCoeffs'
%   Directory containing model coefficients files.
% figDir : char, 1xG, default='figures'
%   Directory to use for output figures.
% figXtn : char, 1xH, default='pdf'
%   Extension to use for figures, which determines the file type.
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% jt_h : double, 1xM, default=[]
%   Ephemeris times of Juno measurements in TDB hours relative to J2000 to use in data comparisons.
% RELATIVE_t : bool, default=0
%   Whether to plot points relative to the start of the input time series. Only has an effect when
%   ``SEQUENTIAL = 0``.
% RELATIVE_r : bool, default=0
%   Whether to plot points relative to distance from the parent planet. Overrides SEQUENTIAL and
%   RELATIVE_t.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 0; end
    if ~exist('coeffpath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('jt_h', 'var'); jt_h = []; end
    if ~exist('RELATIVE_t', 'var'); RELATIVE_t = 0; end
    if ~exist('RELATIVE_r', 'var'); RELATIVE_r = 0; end
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt

    npts = length(t_h);
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    
    if isempty(jt_h); JUNOTOO=0; else; JUNOTOO=1; end
    
    [MagModel, CsheetModel, MPmodel, magModelDescrip, ~] = GetModelOpts(parentName, opt, MPopt);
    magPhase = 0;

    Nmax = 10;
    disp(['Evaluating ' magModelDescrip ' field model with Nmax = ' num2str(Nmax) '.'])
    if contains(magModelDescrip, 'KS2005')
        [Bvec, Mdip_nT, Odip_km] = MagFldJupiterKS2005(r_km, theta, phi, ets, 1);
    else
        [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, r_km, theta, phi, MagModel, ...
            CsheetModel, magPhase, 1, Nmax);
    end
    if ~(strcmp(MPmodel, 'None') || contains(magModelDescrip, 'KS2005'))

        nSW_pcc = 0.14 * ones(1,npts);
        vSW_kms = 400  * ones(1,npts);
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, t_h*3600, xyz_km, Mdip_nT, Odip_km, ...
            parentName, S3coords, MPmodel, coeffPath, 1);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

    end

    Br = Bvec(1,:);
    Bth = Bvec(2,:);
    Bphi = Bvec(3,:);
    BrD = Br - BrSC;
    BthD = Bth - BthSC;
    BphiD = Bphi - BphiSC;
    
    scDataName = [char(scName) ' MAG'];
    commonTitle = [parentName ' model comparison vs ' char(scName) ' data'];
    Rp_km = GetBodyParams(parentName);
    if RELATIVE_r
        [xx, i_xx] = sort(r_km/Rp_km);
        xInfo = [nmTxt 'Distance relative to ' parentName ' center ' mathTxt '(R_P)'];
        fEnd = '_r';
        Br = Br(i_xx);
        BrSC = BrSC(i_xx);
        Bth = Bth(i_xx);
        BthSC = BthSC(i_xx);
        Bphi = Bphi(i_xx);
        BphiSC = BphiSC(i_xx);
        BrD = BrD(i_xx);
        BthD = BthD(i_xx);
        BphiD = BphiD(i_xx);
    else
        fEnd = '';
        if SEQUENTIAL
            xx = 1:npts;
            xInfo = 'Measurement index';
        elseif strcmp(parentName, 'Earth')
            xx = t_h - 175308.0192178;
            xInfo = 'Time relative to NY 2020 (h)';
        elseif strcmp(parentName, 'Uranus')
            xx = t_h + 122154.0036;
            xInfo = 'Time relative to CA (h)';
        elseif strcmp(parentName, 'Neptune')
            xx = t_h + 90752.0566;
            xInfo = 'Time relative to CA (h)';
        elseif RELATIVE_t
            xx = (t_h - t_h(1)) / 24;
            xInfo = 'Time past start (days)';
        else
            xx = t_h;
            xInfo = 'Time past J2000 (h)';
        end
    end

    figNumBase = 3000 + 100*opt + 10*MPopt;
    windowName = [char(scName) 'Br, ' orbStr ', ' magModelDescrip];
    yy = [Br; BrSC];
    yInfo = [mathTxt 'B_r' nmTxt ' component (nT)'];
    legendStrings = [string(magModelDescrip), string(scDataName)];
    titleInfo = commonTitle;
    fName = [char(scName) parentName 'BrComparison' magModelDescrip fEnd];
    PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 1);

    windowName = [char(scName) 'Bth, ' orbStr ', ' magModelDescrip];
    yy = [Bth; BthSC];
    yInfo = [mathTxt 'B_\theta' nmTxt ' component (nT)'];
    fName = [char(scName) parentName 'BthComparison' magModelDescrip fEnd];
    PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 2);

    windowName = [char(scName) 'Bphi, ' orbStr ', ' magModelDescrip];
    yy = [Bphi; BphiSC];
    yInfo = [mathTxt 'B_\phi' nmTxt ' component (nT)'];
    fName = [char(scName) parentName 'BphiComparison' magModelDescrip fEnd];
    PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 3);

    windowName = ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - ' char(scName) ' MAG'];
    yy = [BrD; BthD; BphiD];
    yInfo = 'Component difference (nT)';
    legendStrings = [string([ mathTxt '\Delta B_r']), string([mathTxt '\Delta B_\theta']), ...
        string([mathTxt '\Delta B_\phi'])];
    fName = [char(scName) parentName 'DeltaBComparison' magModelDescrip fEnd];
    PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, fName, ...
        figDir, figXtn, LIVE_PLOTS, figNumBase + 4);

    BrLsq = BrD.^2;
    BthLsq = BthD.^2;
    BphiLsq = BphiD.^2;
    chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3 / npts;
    disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
end
