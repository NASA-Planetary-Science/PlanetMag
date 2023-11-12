% Compare magnetic field measurements from spacecraft against the evaluated
% magnetic field model(s).

% Datasets:
% Swarm MAG data: https://earth.esa.int/web/guest/swarm/data-access
% Current direct link for 1s decimated measurements: https://swarm-diss.eo.esa.int/#swarm/Level1b/Latest_baselines/MAGx_LR

cspice_kclear;
parentName = 'Earth';
scName = "Swarm";
sc = 'Swarm';
SEQUENTIAL = 1; % Whether to plot points by index or hours relative to CA

fullOrbFormatSpec = '%23s%20f%20f%20f%20f%20f%20f%[^\n\r]';
disp(['Importing data files over ' parentName ' orbits.'])
orbStr = [parentName ' data'];
files = dir(fullfile('MAG', sc, '*.tab'));
nTot = 86400*3; % Swarm low-res (LR) decimated data are 1-day files with 1 s cadence, and there are 3 spacecraft
[BrSC, BthSC, BphiSC, r_km, theta, phi] = deal(zeros(1, nTot));
for iFile=1:length(files)
    datFile = files(iFile).name;
    disp(['Loading ' sc ' MAG data from ' datFile '.'])
    fileID = fopen(fullfile('MAG', sc, datFile), 'r');
    magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', 'EndOfLine', '\r\n');
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
    xDescrip = 'Measurement index';
else
    xx = t_h - 175308.0192178;
    xDescrip = 'Time relative to NY 2020 (h)';
end
lat_deg = 90 - rad2deg(theta);
figure; hold on;
set(gcf,'Name', [char(scName) ' latitudes']);
plot(xx, lat_deg, 'DisplayName', 'latitude');
xlabel(xDescrip);
ylabel('Latitude (degrees)');
title('Swarm ABC latitudes on 2020-01-01');
legend();

%% Plot and calculate products
nOpts = 1; nMPopts = 0;
opts = 1:nOpts;
MPopts = -1:-1;
for opt=opts
    for MPopt=MPopts
        GetBplotAndLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
            scName, parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL);
    end
end
