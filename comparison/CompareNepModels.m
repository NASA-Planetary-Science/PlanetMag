% Compare magnetic field measurements from spacecraft against the evaluated
% magnetic field model(s).

% Datasets:
% Voyager 2 MAG: https://doi.org/10.17189/1519975, volume VG2-N-MAG-4-SUMM-NLSCOORDS-12SEC-V1.0

cspice_kclear;
parentName = 'Neptune';
scName = "Voyager 2";
sc = 'Voyager 2';
SEQUENTIAL = 0; % Whether to plot points by index or hours relative to CA

fullOrbFormatSpec = '%23s%10f%8f%8f%10f%10f%10f%[^\n\r]';
disp(['Importing PDS files over ' parentName ' flybys.'])
datFile = fullfile(['MAG/' sc '/vg2_' parentName '_12s_sph.tab']);
orbStr = [parentName ' flyby'];
disp(['Loading ' sc ' MAG data from ' datFile '.'])
fileID = fopen(datFile,'r');
magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', 'EndOfLine', '\r\n');
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
[~, rNep_km] = GetMoonDist(sc, parentName, ets);

% Delete measurement times far from Neptune and junk data
BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
moonProx_RP = 0.1;
PlanetMaxDist_RP = 60;
finiteMax_nT = 17e3;
RPunit = ' R_U';
disp(['Excluding all points satisfying at least one of the following:' newline ...
      'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
      'Suspect measurements, |B| > ' num2str(finiteMax_nT) 'nT.'])
% Full limits
exclude = find(rNep_km/Rp_km > PlanetMaxDist_RP | BmagSC > finiteMax_nT);
ets(exclude) = [];
BrSC(exclude) = [];
BthSC(exclude) = [];
BphiSC(exclude) = [];

npts = length(ets);
t_h = ets / 3600;
disp(['Getting ' sc ' positions for ' num2str(npts) ' pts.'])
[r_km, theta, phi, xyz_km, spkParent] = GetPosSpice(sc, parentName, t_h);


%% Plot and calculate products
nOpts = 1; nMPopts = 1;
opts = 1:nOpts;
MPopts = 1:(nMPopts + 1); % Add 1 to force noMP model in addition
for opt=opts
    for MPopt=MPopts
        GetBplotAndLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
            scName, parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL);
    end
end
