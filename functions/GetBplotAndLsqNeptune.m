function GetBplotAndLsqNeptune(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, scName, ...
    S3coords, orbStr, opt, MPopt, SEQUENTIAL, coeffPath, INC_PDS, INC_O8, INC_PDS_NM3, INC_NM3)
% Plots and calculates comparisons between modeled and measured magnetic fields for Neptune.
%
% Generates time series data of a specified combination of magnetic field models implemented in
% P\lanetMag and compares against spacecraft measurements of the planetary magnetic field.
% Comparisons are plotted and least-squares differences are calculated and printed to the terminal.
% Intended as a final step in model validation; because Voyager 2 trajectories reported in PDS
% files do not perfectly match the trajectories reconstructed from any available SPICE kernels, and
% because the IAU coordinate system for Neptune uses a System II frame (rotating with a stable
% atmospheric feature) instead of a System III frame (rotating with the intrinsic magnetic field),
% this function includes comparisons between evaluation with various relevant coordinate systems.
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
% INC_PDS : bool, default=1
%   Whether to include the full O8 model, including unresolved coefficients, with PDS trajectory in
%   NLS frame in coordinate/frame comparison.
% INC_O8 : bool, default=1
%   Whether to include the full O8 model, including unresolved coefficients, in the NLS frame as
%   implemented in P\lanetMag in coordinate/frame comparison.
% INC_PDS_NM3 : bool, default=1
%   Whether to include O8 model, limited to ``Nmax=3``, with PDS trajectory in NLS frame in
%   coordinate/frame comparison.
% INC_NM3 : bool, default=1
%   Whether to include O8 model, limited to ``Nmax=3``, in NLS frame as implemented in P\lanetMag
%   in coordinate/frame comparison. This implementation uses the Neptune spin pole as referenced at
%   the J2000 epoch and the System III rotation rate as reported in Ness et al. (1989)
%   https://doi.org/10.1126/science.246.4936.1473.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 0; end
    if ~exist('coeffpath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('INC_PDS', 'var'); INC_PDS = 1; end
    if ~exist('INC_O8', 'var'); INC_O8 = 1; end
    if ~exist('INC_PDS_NM3', 'var'); INC_PDS_NM3 = 1; end
    if ~exist('INC_NM3', 'var'); INC_NM3 = 1; end

    parentName = 'Neptune';
    defName = 'O8, NLS relative to J2000 pole (NLS)';
    PDSname = 'O8 with PDS trajectory';
    O8name  = 'O8, NLS with 1989 pole (NLS\_O8)';
    PDS_NM3name = 'O8 with PDS trajectory, Nmax=3';
    NM3name = 'O8, NLS with J2000 pole, Nmax=3';

    npts = length(t_h);
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    
    [MagModel, CsheetModel, MPmodel, magModelDescrip, ~] = GetModelOpts(parentName, opt, MPopt);
    magPhase = 0;

    Nmax = 8;
    disp(['Evaluating ' magModelDescrip ' field model with Nmax = ' num2str(Nmax) '.'])
    if strcmp(magModelDescrip, 'KS2005')
        [Bvec, Mdip_nT, Odip_km] = MagFldJupiterKS2005(r_km, theta, phi, ets, 1);
    else
        [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, r_km, theta, phi, MagModel, ...
            CsheetModel, magPhase, 1, Nmax);
    end
    if ~strcmp(MPmodel, 'None')

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
    
    if INC_PDS || INC_NM3
        BrPDS = Br;
        BthPDS = Bth;
        BphiPDS = Bphi;
        [rNLS, thetaNLS, phiNLS, ~, ~] = GetPosSpice(char(scName), parentName, t_h, 'NLS');
        [BvecNLS, ~, ~] = MagFldParent(parentName, rNLS, thetaNLS, phiNLS, MagModel, ...
            CsheetModel, magPhase, 1, Nmax);
        Br = BvecNLS(1,:);
        Bth = BvecNLS(2,:);
        Bphi = BvecNLS(3,:);
    
        if INC_NM3
            [BvecNM3, ~, ~] = MagFldParent(parentName, rNLS, thetaNLS, phiNLS, MagModel, ...
                CsheetModel, magPhase, 1, 3);
            BrNM3 = BvecNM3(1,:);
            BthNM3 = BvecNM3(2,:);
            BphiNM3 = BvecNM3(3,:);
        end
    end
    
    if INC_PDS_NM3
        [BvecPDS_NM3, ~, ~] = MagFldParent(parentName, r_km, theta, phi, MagModel, CsheetModel, ...
            magPhase, 1, 3);
        BrPDS_NM3 = BvecPDS_NM3(1,:);
        BthPDS_NM3 = BvecPDS_NM3(2,:);
        BphiPDS_NM3 = BvecPDS_NM3(3,:);
    end
    
    if INC_O8
        [rO8, thetaO8, phiO8, ~, ~] = GetPosSpice(char(scName), parentName, t_h, 'NLS_RADEC');
        [BvecO8, ~, ~] = MagFldParent(parentName, rO8, thetaO8, phiO8, MagModel, CsheetModel, ...
            magPhase, 1, Nmax);
        BrO8 = BvecO8(1,:);
        BthO8 = BvecO8(2,:);
        BphiO8 = BvecO8(3,:);
    end
    
    commonTitle = ['Neptune field model comparison, Nmax=' num2str(Nmax)];
    if SEQUENTIAL
        xx = 1:npts;
        xDescrip = 'Measurement index';
    else
        xx = t_h + 90752.0566;
        xDescrip = 'Time relative to CA (h)';
    end
    
    figure; hold on;
    set(gcf,'Name', ['Br, ' orbStr ', ' magModelDescrip]);
    plot(xx, Br, 'DisplayName', defName);
    if INC_PDS; plot(xx, BrPDS, 'DisplayName', PDSname); end
    if INC_NM3; plot(xx, BrNM3, 'DisplayName', NM3name); end
    if INC_PDS_NM3; plot(xx, BrPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_O8;  plot(xx, BrO8, 'DisplayName', O8name); end
    plot(xx, BrSC, 'DisplayName', 'Voyager 2 MAG');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    xlim([-1.5,1.5])
    legend();

    figure; hold on;
    set(gcf,'Name', ['Bth, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bth, 'DisplayName', defName);
    if INC_PDS; plot(xx, BthPDS, 'DisplayName', PDSname); end
    if INC_NM3; plot(xx, BthNM3, 'DisplayName', NM3name); end
    if INC_PDS_NM3; plot(xx, BthPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_O8;  plot(xx, BthO8, 'DisplayName', O8name); end
    plot(xx, BthSC, 'DisplayName', 'Voyager 2 MAG');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    xlim([-1.5,1.5])
    legend();

    figure; hold on;
    set(gcf,'Name', ['Bphi, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bphi, 'DisplayName', defName);
    if INC_PDS; plot(xx, BphiPDS, 'DisplayName', PDSname); end
    if INC_PDS_NM3; plot(xx, BphiPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_NM3; plot(xx, BphiNM3, 'DisplayName', NM3name); end
    if INC_O8;  plot(xx, BphiO8, 'DisplayName', O8name); end
    plot(xx, BphiSC, 'DisplayName', 'Voyager 2 MAG');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    xlim([-1.5,1.5])
    legend();
    
    % Evaluate magnitudes
    Bmag = sqrt(Br.^2 + Bth.^2 + Bphi.^2);
    BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
    if INC_PDS; BmagPDS = sqrt(BrPDS.^2 + BthPDS.^2 + BphiPDS.^2); end
    if INC_PDS_NM3; BmagPDS_NM3 = sqrt(BrPDS_NM3.^2 + BthPDS_NM3.^2 + BphiPDS_NM3.^2); end
    if INC_NM3; BmagNM3 = sqrt(BrNM3.^2 + BthNM3.^2 + BphiNM3.^2); end
    if INC_O8;  BmagO8 = sqrt(BrO8.^2 + BthO8.^2 + BphiO8.^2); end
    
    % Plot with same axes as past comparisons
    figure; hold on;
    set(gcf,'Name', ['Bmag, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bmag, 'DisplayName', defName);
    if INC_PDS; plot(xx, BmagPDS, 'DisplayName', PDSname); end
    if INC_PDS_NM3; plot(xx, BmagPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_NM3; plot(xx, BmagNM3, 'DisplayName', NM3name); end
    if INC_O8;  plot(xx, BmagO8, 'DisplayName', O8name); end
    plot(xx, BmagSC, 'DisplayName', 'Voyager 2 MAG');
    set(gca, 'YScale', 'log');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    legend();
    grid on;
    xlim([-1.5,1.5])
    ylim([1e2,1e5])
    
    % Plot with same axes as Connerney et al. (1991)
    figure; hold on;
    set(gcf,'Name', ['Bmag (C1991), ' orbStr ', ' magModelDescrip]);
    plot(xx, Bmag, 'DisplayName', defName);
    if INC_PDS; plot(xx, BmagPDS, 'DisplayName', PDSname); end
    if INC_PDS_NM3; plot(xx, BmagPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_NM3; plot(xx, BmagNM3, 'DisplayName', NM3name); end
    if INC_O8;  plot(xx, BmagO8, 'DisplayName', O8name); end
    plot(xx, BmagSC, 'DisplayName', 'Voyager 2 MAG');
    set(gca, 'YScale', 'log');
    xlabel('Time of day 1989 Aug 25');
    delt = 1 - (55*60+40.076)/3600;
    xticks([(-6:-1)/6 + delt, 0, (0:6)/6 + delt])
    xticklabels({'03:00','03:10','03:20','03:30','03:40','03:50','CA','04:00','04:10','04:20', ...
        '04:30','04:40','04:50','05:00'})
    ylabel('Vector component (nT)');
    title('Field magnitude comparisons, axes as in C1991');
    legend();
    grid on;
    xlim([-1+delt 1+delt])
    ylim([1e2,1e5])

    BrD = Br - BrSC;
    BthD = Bth - BthSC;
    BphiD = Bphi - BphiSC;
    BmagD = Bmag - BmagSC;

    figure; hold on;
    set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - MAG, ' defName]);
    plot(xx, BrD);
    plot(xx, BthD);
    plot(xx, BphiD);
    plot(xx, BmagD);
    xlabel(xDescrip);
    ylabel('Component difference (nT)');
    xlim([-1,1])
    title('Neptune field model comparison for Nmax=8, NLS relative to J2000 pole - Voyager 2 MAG')
    legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|')
    
    BrLsq = BrD.^2;
    BthLsq = BthD.^2;
    BphiLsq = BphiD.^2;
    chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3/npts;
    disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
    
    %%
    if INC_PDS
        BrD = BrPDS - BrSC;
        BthD = BthPDS - BthSC;
        BphiD = BphiPDS - BphiSC;
        BmagD = BmagPDS - BmagSC;

        figure; hold on;
        set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ', PDS trajec - MAG']);
        plot(xx, BrD);
        plot(xx, BthD);
        plot(xx, BphiD);
        plot(xx, BmagD);
        xlabel(xDescrip);
        ylabel('Component difference (nT)');
        xlim([-1,1])
        title('Neptune field model comparison, PDS trajec - MAG')
        legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|')
        
        BrLsq = BrD.^2;
        BthLsq = BthD.^2;
        BphiLsq = BphiD.^2;
        chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3/npts;
        disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
    end
    
    %%
    if INC_O8
        BrD = BrO8 - BrSC;
        BthD = BthO8 - BthSC;
        BphiD = BphiO8 - BphiSC;
        BmagD = BmagO8 - BmagSC;

        figure; hold on;
        set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - MAG, ' O8name]);
        plot(xx, BrD);
        plot(xx, BthD);
        plot(xx, BphiD);
        plot(xx, BmagD);
        xlabel(xDescrip);
        ylabel('Component difference (nT)');
        xlim([-1,1])
        title('Neptune field model comparison, NLS exactly as defined with O8 model');
        legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|');
        
        BrLsq = BrD.^2;
        BthLsq = BthD.^2;
        BphiLsq = BphiD.^2;
        chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3/npts;
        disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
    end
end
