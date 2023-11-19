function CompareUraModels(LIVE_PLOTS, scName, SEQUENTIAL)
% Compare magnetic field measurements from spacecraft near Uranus against each implemented magnetic
% field model.
%
% Datasets:
%
%   * Voyager 2 MAG: https://doi.org/10.17189/1520034, volume VG2-U-MAG-4-RDR-U1COORDS-48SEC-V1.0
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% scName : string, default="Voyager 2"
%   Spacecraft name for which magnetic field data will be compared against implemented models. A
%   directory must exist with this name in the ``MAG`` directory within the top-level P\lanetMag
%   directory. This directory will be searched for body-specific data files, and each of these
%   files will be loaded.
% SEQUENTIAL : bool, default=0
%   Whether to plot points by index or hours relative to a reference time (typically closest
%   approach).

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('scName', 'var'); scName = "Voyager 2"; end
    if ~exist('SEQUENTIAL', 'var'); SEQUENTIAL = 0; end
    
    cspice_kclear;
    parentName = 'Uranus';
    sc = char(scName);
    
    fullOrbFormatSpec = '%24s%5d%2d%2d%1d%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%3d%8s%[^\n\r]';
    disp(['Importing PDS files over ' parentName ' flybys.'])
    datFile = fullfile(['MAG/' sc '/vg2_' parentName '_48s_sph.tab']);
    orbStr = [parentName ' flyby'];
    disp(['Loading ' sc ' MAG data from ' datFile '.'])
    fileID = fopen(datFile,'r');
    magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', ',', 'TextType', 'char', 'EndOfLine', '\r\n');
    fclose(fileID);
    
    t_UTC = magData{1}';
    BrSC = magData{6}';
    BthSC = magData{7}';
    BphiSC = magData{8}';
    
    LoadSpice(parentName, sc);
    Rp_km = 25559;    
    
    nTot = length(t_UTC);
    disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
    ets = cspice_str2et(t_UTC);
    disp(['Getting moon and planet distances for all ' num2str(nTot) ' points.'])
    [~, rP_km] = GetMinMoonDist(sc, parentName, ets);
    
    % Delete measurement times far from Uranus and junk data
    BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
    moonProx_RP = 0.1;
    PlanetMaxDist_RP = 60;
    finiteMax_nT = 17e3;
    RPunit = [' R_' parentName(1)];
    disp(['Excluding all points satisfying at least one of the following:' newline ...
          'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
          'Suspect measurements, |B| > ' num2str(finiteMax_nT) 'nT.'])
    % Full limits
    exclude = find(rP_km/Rp_km > PlanetMaxDist_RP | BmagSC > finiteMax_nT);
    ets(exclude) = [];
    BrSC(exclude) = [];
    BthSC(exclude) = [];
    BphiSC(exclude) = [];
    
    npts = length(ets);
    t_h = ets / 3600;
    disp(['Getting ' sc ' positions for ' num2str(npts) ' pts.'])
    [r_km, theta, phi, xyz_km, spkParent] = GetPosSpice(sc, parentName, t_h);
    
    
    %% Plot and calculate products
    nOpts = 2; nMPopts = 2;
    opts = 1:nOpts;
    MPopts = 1:(nMPopts + 1); % Add 1 to force noMP model in addition
    for opt=opts
        for MPopt=MPopts(2:end)
            GetBplotAndLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
                scName, parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL);
        end
    end
        
    
    %% Plot trajectory
    t_h = linspace(cspice_str2et('1986-01-23T00:00:00.000'), cspice_str2et('1986-01-27T00:00:00.000'), 5000)/3600;
    [~, ~, ~, xyzUSO_km, ~] = GetPosSpice(sc, parentName, t_h, 'USO');
    xyz_Rp = xyzUSO_km / Rp_km;
    x = xyz_Rp(1,:); y = xyz_Rp(2,:); z = xyz_Rp(3,:);
    clf(); hold on;
    plot3(x,y,z, 'LineWidth', 1.5);
    the = linspace(0,pi,50); ph = linspace(0,2*pi,100);
    [the2D, ph2D] = meshgrid(the,ph);
    xp = sin(the2D) .* cos(ph2D); yp = sin(the2D) .* sin(ph2D); zp = cos(the2D);
    surf(xp,yp,zp, 'FaceColor', 'c');
    axlim = 30;
    pbaspect([1 1 1]);
    xlim([-axlim,axlim]);ylim([-axlim,axlim]);zlim([-axlim,axlim]);
    grid on;
    xlabel('x USO (R_U, toward Sun)');ylabel('y USO (R_U), \approx orbital v');zlabel('z USO (R_U)');
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Ydir', 'reverse')
    
    etN = cspice_str2et('1986-01-24T17:58:51.346');
    rotMat = cspice_pxform('IAU_URANUS', 'USO', etN);
    vecMat = zeros(3,1);
    vecMat(3,1) = -2.5;
    poleUSO = squeeze(pagemtimes(rotMat, vecMat));
    plot3([0,poleUSO(1)], [0,poleUSO(2)], [0,poleUSO(3)], 'Color', 'r', 'LineWidth', 2)
    scatter3(poleUSO(1), poleUSO(2), poleUSO(3), 15, 'r')
    
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
    
    
    %% Plot L shell
    t_h = linspace(cspice_str2et('1986-01-23T00:00:00.000'), cspice_str2et('1986-01-27T00:00:00.000'), 500)/3600;
    [r_km, theta, ~, ~, ~] = GetPosSpice(sc, parentName, t_h, 'SMU');
    lambda = pi - theta;
    r_Rp = r_km / Rp_km;
    L = r_Rp ./ cos(lambda).^2;
    upperLim = 40000;
    DO_MPAUSE = 0;
    
    figure (100); clf(); hold on;
    plot(t_h + 122154.0036, L, 'DisplayName', 'L shell', 'LineWidth', 2)
    xlabel('Time relative to Uranus CA (h)');
    ylabel('L shell');
    ylim([0 upperLim])
    set(gca, 'YScale', 'log')
    title('Uranus L shell during Voyager 2 flyby');
    
    figure (101); clf(); hold on;
    plot(r_Rp, L, 'DisplayName', 'L shell', 'LineWidth', 2)
    xlabel('r (R_U)');
    ylabel('L shell');
    ylim([0 upperLim])
    set(gca, 'YScale', 'log')
    title('Uranus L shell during Voyager 2 flyby vs. radial distance');

end