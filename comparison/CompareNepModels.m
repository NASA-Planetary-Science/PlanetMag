function CompareNepModels(LIVE_PLOTS, scName, SEQUENTIAL)
% Compare magnetic field measurements from spacecraft near Neptune against each implemented
% magnetic field model.
%
% Datasets:
%
%   * Voyager 2 MAG: https://doi.org/10.17189/1519975, volume VG2-N-MAG-4-SUMM-NLSCOORDS-12SEC-V1.0
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
    parentName = 'Neptune';
    sc = char(scName);
    
    fullOrbFormatSpec = '%23s%10f%8f%8f%10f%10f%10f%[^\n\r]';
    disp(['Importing PDS files over ' parentName ' flybys.'])
    datFile = fullfile(['MAG/' sc '/vg2_' parentName '_12s_sph.tab']);
    orbStr = [parentName ' flyby'];
    disp(['Loading ' sc ' MAG data from ' datFile '.'])
    fileID = fopen(datFile,'r');
    magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', ...
        'EndOfLine', '\r\n');
    fclose(fileID);
    
    t_UTC = magData{1}';
    BrSC = magData{5}';
    BthSC = magData{6}';
    BphiSC = magData{7}';
    
    LoadSpice(parentName, sc);
    Rp_km = 24765;    
    
    nTot = length(t_UTC);
    disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
    ets = cspice_str2et(t_UTC);
    disp(['Getting moon and planet distances for all ' num2str(nTot) ' points.'])
    [~, rP_km] = GetMinMoonDist(sc, parentName, ets);
    
    % Delete measurement times far from Neptune and junk data
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
    
    tCA_h = -90752.0566;
    tRel_h = t_h - tCA_h;
    Descrip = 'Time relative to CA (h)';
    
    COMPARE_PDS = 1;
    if COMPARE_PDS
        offset = 0;%-0.0825;
        offset_th = 0;
        r = magData{2}' * Rp_km; r(exclude) = [];
        th = deg2rad(offset_th + 90 - magData{3}'); th(exclude) = [];
        ph = deg2rad(offset -magData{4}'); ph(exclude) = [];
        r_Rp = r / Rp_km;
        x = r .* sin(th) .* cos(ph);
        y = r .* sin(th) .* sin(ph);
        z = r .* cos(th);
        figure; hold on
        dx = x - xyz_km(1,:);
        dy = y - xyz_km(2,:);
        dz = z - xyz_km(3,:);
        dr = sqrt(dx.^2 + dy.^2 + dz.^2);
        plot(tRel_h, dx, 'DisplayName', '\Delta x')
        plot(tRel_h, dy, 'DisplayName', '\Delta y')
        plot(tRel_h, dz, 'DisplayName', '\Delta z')
        plot(tRel_h, dr, 'DisplayName', '\Delta r')
        xlabel(Descrip);
        ylabel('Coordinate diff (km)');
        title('Location difference in NLS frame, PDS - SPICE');
        legend()
        r_km = r; theta = th; phi = ph; xyz_km = [x; y; z];
    end
    
    
    %% Plot and calculate products
    nOpts = 1; nMPopts = 0;
    opts = 1:nOpts;
    %MPopts = 1:(nMPopts + 1); % Add 1 to force noMP model in addition
    MPopts = -1:-1;
    for opt=opts
        for MPopt=MPopts
            GetBplotAndLsqNeptune(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
                scName, spkParent, orbStr, opt, MPopt, SEQUENTIAL, 1, 1, 1, 1);
        end
    end
    
    
    %% Plot trajectory
    figure(1000); clf(); hold on;
    set(gcf,'Name', 'Trajectories');
    the = linspace(0,pi,19); ph = linspace(0,2*pi,37);
    [the2D, ph2D] = meshgrid(the,ph);
    xp = sin(the2D) .* cos(ph2D); yp = sin(the2D) .* sin(ph2D); zp = cos(the2D);
    surf(xp,yp,zp, 'FaceColor', 'b');
    plot3([0,1.25], [0,0], [0,0], 'Color', 'r', 'LineWidth', 2)
    scatter3(1.25, 0, 0, 15, 'r')
    axlim = 1.5;
    pbaspect([1 1 1]);
    xlim([-axlim,axlim]);ylim([-axlim,axlim]);zlim([-axlim,axlim]);
    grid on;
    xlabel('Despun NLS x (R_S)');ylabel('Despun NLS y (R_S)');zlabel('Despun NLS z (R_S)');
    title('Voyager 2 Neptune trajectory in despun NLS')
    
    tNLS_h = -90752.05106;
    tOff_h = tNLS_h - tCA_h;
    despin = 2*pi * (tRel_h - tOff_h) / 16.11;
    xSC = r_Rp .* sin(theta) .* cos(phi + despin);
    ySC = r_Rp .* sin(theta) .* sin(phi + despin);
    zSC = r_Rp .* cos(theta);
    scTraj{2} = plot3(xSC,ySC,zSC, 'LineWidth', 1.5); name{2} = "PDS";
    rCA = 1.180899;
    thCA = deg2rad(12.45913);
    phiCA = deg2rad(-167.7);
    CA = [rCA*sin(thCA)*cos(phiCA), rCA*sin(thCA)*sin(phiCA), rCA*cos(thCA)];
    CAlineSC = plot3([0,CA(1)], [0,CA(2)], [0,CA(3)], 'Color', 'k', 'LineWidth', 2);
    CAnameSC = "CA";
    scatter3(CA(1), CA(1), CA(1), 15, 'k')
    
    [~, ~, ~, xyz_km, ~] = GetPosSpice(sc, parentName, t_h, 'NLS');
    xyz_Rp = xyz_km / Rp_km;
    r_Rp = sqrt(xyz_Rp(1,:).^2 + xyz_Rp(2,:).^2 + xyz_Rp(3,:).^2);
    theta = acos(xyz_Rp(3,:) ./ r_Rp);
    phi = atan2(xyz_Rp(2,:), xyz_Rp(1,:));
    
    x = r_Rp .* sin(theta) .* cos(phi + despin);
    y = r_Rp .* sin(theta) .* sin(phi + despin);
    z = r_Rp .* cos(theta);
    scTraj{1} = plot3(x,y,z, 'LineWidth', 1.5); name{1} = "SPICE in NLS";
    
    [~, ~, ~, CAN, ~] = GetPosSpice(sc, parentName, tNLS_h, 'NLS');
    CAN = CAN / Rp_km;
    CAline = plot3([0,CAN(1)], [0,CAN(2)], [0,CAN(3)], 'Color', 'y', 'LineWidth', 2);
    CAname = "CA (NLS)";
    scatter3(CAN(1), CAN(1), CAN(1), 15, 'y')
    
    [~, ~, ~, xyz_km, ~] = GetPosSpice(sc, parentName, t_h, 'NLS_RADEC');
    [x, y, z] = GetDespun(xyz_km/Rp_km, despin);
    scTraj{3} = plot3(x,y,z, 'LineWidth', 1.5, 'Color', 'm');
    name{3} = "SPICE in NLS as defined in O8";
    
    legend([scTraj{1} scTraj{2} scTraj{3}], [name{1}, name{2}, name{3}])

end
