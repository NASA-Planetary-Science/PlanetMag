cspice_kclear;
moonName = 'Ganymede';
parentName = 'Jupiter';
scName = "Galileo";
orbNum = -2;
SEQUENTIAL = 1; % Whether to plot points by index or hours past J2000
FULLORBITS = 1; % Whether to evalate goodness of fit from full-orbit data in System III coordinates
FLYBYS = 1; % Whether to evaluate goodness of fit from flyby data in moon coordinates
JUNOTOO = 1;

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

if FULLORBITS
    [t_UTC, BrSC, BthSC, BphiSC] = deal([]);
    fullOrbFormatSpec = '%23s%24s%10f%10f%10f%10f%7f%7f%7f%f%[^\n\r]';
    disp(['Importing PDS files over ' orbStr '.'])
    for i=1:length(orbList)
        if orbNum ~= -1
            disp(['Loading ' char(scName) ' MAG data for orbit ' num2str(orbList(i)) '.'])
        end
        datFile = fullfile(['MAG/' char(scName) '/Jupiter/ORB' sprintf('%02d', orbList(i)) '_SYS3.TAB']);
        fileID = fopen(datFile,'r');
        magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', 'EndOfLine', '\r\n');
        fclose(fileID);

        t_UTC = [t_UTC, magData{1}'];
        BrSC = [BrSC, magData{3}'];
        BthSC = [BthSC, magData{4}'];
        BphiSC = [BphiSC, magData{5}'];
    end

    LoadSpice(moonName, char(scName));
    sc = 'Galileo Orbiter';
    [R_P, ~, ~, ~, ~, ~, ~] = GetBodyParams(moonName);

    nTot = length(t_UTC);
    disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
    ets = cspice_str2et(t_UTC);
    disp(['Getting moon- and Jupiter-relative distances for all ' num2str(nTot) ' points.'])
    [rMinMoon_km, rJup_km] = GetMoonDist(sc, parentName, ets);

    % Delete measurement times near moons, far from Jupiter, and junk data
    BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
    moonProx_RP = 0.1;
    PlanetMinDist_RP = 10;
    PlanetMaxDist_RP = 20;
    finiteMax_nT = 2500;
    RPunit = ' R_J';
    disp(['Excluding all points satisfying at least one of the following:' newline ...
          'Distance to a major moon < ' num2str(moonProx_RP) RPunit newline ...
          'Planetocentric distance < ' num2str(PlanetMinDist_RP) RPunit newline ...
          'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
          'Suspect measurements, |B| > ' num2str(finiteMax_nT) 'nT.'])
    % Full limits
    exclude = find(rMinMoon_km/R_P < moonProx_RP | rJup_km/R_P < PlanetMinDist_RP | rJup_km/R_P > PlanetMaxDist_RP | BmagSC > finiteMax_nT);
    % No parent planet distance limits (full magnetosphere)
    %exclude = find(rMinMoon_km/R_P < moonProx_RP | BmagSC > finiteMax_nT);
    ets(exclude) = [];
    BrSC(exclude) = [];
    BthSC(exclude) = [];
    BphiSC(exclude) = [];

    npts = length(ets);
    t_h = ets / 3600;
    disp(['Getting ' sc ' positions for ' num2str(npts) ' pts over ' orbStr '.'])
    [r_km, lat_deg, lon_deg] = GetPosSpice(sc, parentName, t_h);
    alt_km = r_km - R_P;
    theta = deg2rad(90 - lat_deg);
    phi = deg2rad(lon_deg);
end
spkParent = ['IAU_' upper(parentName)];
spkMoon = ['IAU_' upper(moonName)];

%% Flyby data

if FLYBYS
    [fbt_UTC, fbBrSC, fbBthSC, fbBphiSC] = deal([]);
    flybyFormatSpec = '%23s%10f%10f%10f%10f%7f%7f%7f%f%[^\n\r]';
    for i=1:length(fbList)
        if orbNum ~= -1
            disp(['Loading ' char(scName) ' MAG data for flyby ' moonName(1) num2str(orbList(i)) '.'])
        end
        datFile = fullfile(['MAG/' char(scName) '/' moonName '/ORB' sprintf('%02d_', fbList(i)) fbCode '_SYS3.TAB']);
        fileID = fopen(datFile,'r');
        magData = textscan(fileID, flybyFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', 'EndOfLine', '\r\n');
        fclose(fileID);

        fbt_UTC = [fbt_UTC, magData{1}'];
        fbBrSC = [fbBrSC, magData{2}'];
        fbBthSC = [fbBthSC, magData{3}'];
        fbBphiSC = [fbBphiSC, magData{4}'];
    end

    fbets = cspice_str2et(fbt_UTC);
    fbt_h = fbets / 3600;
    [fbr_km, fblat_deg, fblon_deg] = GetPosSpice(sc, parentName, fbt_h);
    fbalt_km = fbr_km - R_P;
    fbtheta = deg2rad(90 - fblat_deg);
    fbphi = deg2rad(fblon_deg);
    
    [BxSCS3, BySCS3, BzSCS3] = Bsph2Bxyz(fbBrSC, fbBthSC, fbBphiSC, fbtheta, fbphi);
    [BxSC, BySC, BzSC] = RotateBspice(BxSCS3, BySCS3, BzSCS3, fbets, spkParent, spkMoon);
    r_RM = GetTargetMoonDist(sc, moonName, parentName, fbets);
end

% Add Juno flyby data
if JUNOTOO && FLYBYS
    sc = 'JUNO';
    cspice_kclear;
    LoadSpice(moonName, 'Juno');
    fbjNum = 34;
    datFile = fullfile(['MAG/Juno/' moonName '/ORB' num2str(fbjNum) '_FGM_IAU.txt']);
    disp(['Loading FGM data for Juno ' moonName ' flyby on orbit ' num2str(fbjNum) '.'])
    junoFlybyFormatSpec = '%23s%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(datFile,'r');
    magData = textscan(fileID, junoFlybyFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', 'EndOfLine', '\r\n');
    
    jt_UTC = magData{1}';
    jBx = magData{2}';
    jBy = magData{3}';
    jBz = magData{4}';
    jets = cspice_str2et(jt_UTC);
    jt_h = jets / 3600;
    
    [jr_km, jlat_deg, jlon_deg] = GetPosSpice(sc, parentName, jt_h);
    jalt_km = jr_km - R_P;
    jtheta = deg2rad(90 - jlat_deg);
    jphi = deg2rad(jlon_deg);
    jr_RM = GetTargetMoonDist(sc, moonName, parentName, jets);
    
    % Append to lists of coords to eval in Jupiter field model
    fbets = [fbets, jets];
    % NOT t_h, we will use the separate arrays to organize outputs.
    fbalt_km = [fbalt_km, jalt_km];
    fblat_deg = [fblat_deg, jlat_deg];
    fblon_deg = [fblon_deg, jlon_deg];
    r_RM = [r_RM, jr_RM];
    BxSC = [BxSC, jBx];
    BySC = [BySC, jBy];
    BzSC = [BzSC, jBz];
    
    scName = [scName, "Juno"];
else
    jt_h = [];
end

%% Plot and calculate products

for opt=5:6
    if FULLORBITS
        GetBplotAndLsq(ets, t_h, alt_km, lat_deg, lon_deg, BrSC, BthSC, BphiSC, ...
            scName, parentName, orbStr, opt, SEQUENTIAL);
    end
    
    if FLYBYS
        GetBplotAndLsqMoon(fbets, fbt_h, fbalt_km, fblat_deg, fblon_deg, ...
            r_RM, BxSC, BySC, BzSC, ...
            scName, parentName, moonName, era, fbStr, opt, SEQUENTIAL, jt_h);
    end
end
