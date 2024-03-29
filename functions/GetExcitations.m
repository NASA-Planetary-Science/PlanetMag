function fList = GetExcitations(moonName, etMid_day)
% Retrieve a list of excitation frequencies to invert.
%
% Returned excitation frequencies consist of combinations of precise orbital periods, including
% those adjusted for nutation and precession. For most bodies, one or more empirically determined
% oscillation periods are also included, such as those corresponding to the true anomaly period.
% Only one of the parameters marked "adjusted" should be used for each body.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of moon for which to return excitation frequencies.
% etMid_day : double, default=0
%   Ephemeris time (ET) near the center of time series to be evaluated in days. Default is J2000.
%
% Returns
% -------
% fList : table, P'x2
%   Contains columns:
%   
%   - ``double``: Frequencies of excitations in Hz.
%   - ``string``: Descriptive names for each excitation.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('etMid_day', 'var'); etMid_day = 0; end

    fEyr = 1 / 365.256 / 86400; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    fJyr = 1 / 4332.589 / 86400; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
    fSyr = 1 / 10759.22 / 86400; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
    fUyr = 1 / 30685.4 / 86400; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranusfact.html
    fNyr = 1 / 60189 / 86400; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html
    % Fetch remaining quantities from PCK file
    [~, ~, ~, ~, omegaParent_radps, omegaMoon_radps, ~, ~, nutPrecParent, nutPrecMoon] ...
        = GetBodyParams(moonName);
    fPar = abs(omegaParent_radps / 2/pi);
    fOrb = abs(omegaMoon_radps / 2/pi);
    Pnut0 = nutPrecParent(1:2:end);
    Pnut1 = nutPrecParent(2:2:end) / 36525;
    fOrbAdj = abs(omegaMoon_radps / 2/pi + sum(nutPrecMoon .* Pnut1 .* cosd(Pnut0 ...
        + Pnut1*etMid_day)) / 360 / 86400);
    fSyn = fPar - fOrb;
    fSynAdj = fPar - fOrbAdj;
    
    switch moonName

        %% Earth moon
        case 'Moon'

            fList = [
                {fSyn, "synodic"}
                {fOrb, "orbital"}
                {2*fSyn, "synodic 2nd"}
                ];

        %% Jupiter moons
        case 'Io'
            fTA = 1 / 3600 / 42.30632; % True anomaly period

            fList = [
                {fSyn, "synodic"}
                {fTA, "true anomaly"}
                {2*fSyn, "synodic 2nd"}
                {3*fSyn, "synodic 3rd"}
                {4*fSyn, "synodic 4th"}
                {fOrbAdj + fJyr, "orbital+year beat"}
                {fSynAdj - fTA, "synodic-TA beat"}
                {fSynAdj + fTA, "synodic+TA beat"}
                {2*fSynAdj - fTA, "synodic 2nd-TA beat"}
                ];

        case 'Europa'
            fTA = 1 / 84.63681 / 3600;
             
            fList = [
                {fSyn, "synodic"}
                {fTA, "true anomaly"}
                {fOrbAdj, "orbital"}
                {2*fSyn, "synodic 2nd"}
                {3*fSyn, "synodic 3rd"}
                {fSynAdj - fTA, "synodic-TA beat"}
                {fSynAdj + fTA, "synodic+TA beat"}
                {2*fSyn - fTA, "synodic 2nd-TA beat"}
                {fTA - fJyr, "TA-year beat"}
                {fTA + fJyr, "TA+year beat"}
                {fOrbAdj - fJyr, "orbital-year beat"}
                ];

        case 'Ganymede'
            fList = [
                {fSyn,    "synodic"}
                {fOrbAdj, "orbital"}
                {2*fSyn,  "synodic 2nd"}
                ];

        case 'Callisto'

            fList = [
                {fSyn, "synodic"}
                {fOrb, "orbital"}
                {2*fSyn, "synodic 2nd"}
                {3*fSyn, "synodic 3rd"}
                {5*fSyn, "synodic 5th"}
                ];

        %% Saturn moons         
        case 'Mimas'
            fTA = 1 / 3600 / 22.55972;
            fList = [
                {fTA, "true anomaly"}
                {fOrb, "orbital"}
                {2*fOrb, "orbital 2nd"}
                {2*fOrb - fTA, "orbital 2nd-TA beat"}
                {2*fTA - fOrb, "TA 2nd-orbital beat"}
                {2*fOrb - fTA - fSyr, "orbital 2nd-TA-year beat"}
                {2*fOrb - fTA - 2*fSyr, "orbital 2nd-TA-half year beat"}
                {4*fOrb - 2*fTA - 2*fSyr, "orbital 4th-TA 2nd-half year beat"}
                ];

        case 'Enceladus'
            fTA = 1 / 3600 / 32.92759; % True anomaly period
            fList = {fTA, "true anomaly"};

        case 'Tethys'
            fTA = 1 / 3600 / 45.25957; % True anomaly period
            fList = [
                {fTA, "true anomaly"}
                {fTA - fSyr, "TA-year beat"}
                ];

        case 'Dione'            
            fList = {fOrbAdj - 2*fSyr, "orbital-half year beat"};

        case 'Rhea'
            fList = {fOrb, "orbital"};

        case 'Titan'
            fList = {fOrb, "orbital"};

        case 'Iapetus'
            fList = {fOrb, "orbital"};

        %% Uranus moons   
        case 'Miranda'
            fList = [
                {fSyn, "synodic"}
                {fOrb, "orbital"}
                {2*fSyn, "synodic 2nd"}
                {3*fSyn, "synodic 3rd"}
                {4*fSyn, "synodic 4th"}
                {fOrb - fSyn, "orbital-synodic beat"}
                {fSyn + fOrb, "synodic+orbital beat"}
                {2*fSyn + fOrb, "synodic 2nd+orbital beat"}
                ];

        case 'Ariel'
            fList = [
                {fSyn, "synodic"}
                {2*fSyn, "synodic 2nd"}
                {3*fSyn, "synodic 3rd"}
                ];

        case 'Umbriel'
            fList = [
                {fSyn, "synodic"}
                {2*fSyn, "synodic 2nd"}
                {fSyn - fOrb, "synodic-orbital beat"}
                ];

        case 'Titania'
            fList = [
                {fSyn, "synodic"}
                {2*fSyn, "synodic 2nd"}
                ];

        case 'Oberon'
            fList = {fSyn, "synodic"};


        %% Neptune moons
        case 'Triton'
            fSyn = fPar + fOrb; % Retrograde orbit
            fList = [
                {fSyn, "synodic"}
                {fOrb, "orbital"}
                {fSyn - fOrb, "synodic-orbital beat"}
                ];
                    
    end
    
    fList = sortrows(cell2table(fList), 1, 'descend');
    if any(contains(fList{:,2}, ','))
        error(['At least one excitation period name contains a comma for ' moonName ...
            ':' newline sprintf('\n%s', fList{:,2}) '\nRemove the comma to be compatible ' ...
            'with .csv file output.'])
    end
end
