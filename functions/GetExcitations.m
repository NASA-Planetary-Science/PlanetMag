function [f, fNames] = GetExcitations(moonName, etMid_day)
% Retrieve a list of excitation frequencies to invert.
%
% Returned excitation frequencies consist of combinations of precise orbital periods, including
% those adjusted for nutation and precession. For most bodies, one or more empirically determined
% oscillation periods are also included, such as those corresponding to the true anomaly period.
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
% f : double, 1xP'
%   Frequencies of excitations in Hz.
% fNames : string, 1xP'
%   Descriptive names for each excitation.

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

            f = [fSyn, fOrb];

            f = [f, 2*fSyn, ...
                ];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                ];

        %% Jupiter moons 
        case 'Io'
            fIoTA = 1 / 3600 / 42.315044531808887029; % true anomaly period
            fOrbMid1 = 1 / 3600 / 42.431381950208163;
            fOrbMid2 = 1 / 3600 / 42.305626942856122;

            f = [fSyn, fOrbMid1, fOrbMid2];

            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - fIoTA, ... % 1st harmonic beats
                    fSyn + fIoTA, ...
                    2*fSyn - fIoTA, ... % 2nd harmonic beats
                    fSynAdj - fIoTA, ... % 1st harmonic beats
                    fSynAdj + fIoTA, ...
                    2*fSynAdj - fIoTA, ... % 2nd harmonic beats
                    fSynAdj - fOrbMid2, ...  % Beat between strong(est) oscillations
                    ];

            fNames = [
                "synodic"
                "orbital mid 1"
                "orbital mid 2"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic-TA beat"
                "synodic+TA beat"
                "synodic 2nd-TA beat"
                "adjusted synodic-TA beat"
                "adjusted synodic+TA beat"
                "adjusted synodic 2nd-TA beat"
                "adjusted synodic-orbital mid 2 beat"
                ];

        case 'Europa'
            fTA = 1 / 84.62749 / 3600;
             
            f = [fSyn, fTA, fOrbAdj];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    fSyn - fTA, ... % 1st harmonic beats
                    fSyn + fTA, ...
                    fTA - fJyr, ... % Solar beats
                    fTA + fJyr, ... % with orbital
                    fOrbAdj - fJyr, ... % periods
                    ];

            fNames = [
                "synodic"
                "true anomaly"
                "adjusted orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic-TA beat"
                "synodic+TA beat"
                "TA-year beat"
                "TA+year beat"
                "adjusted orbital-year beat"
                ];

        case 'Ganymede'
            f = [fSyn, fOrbAdj];

            f = [f, 2*fSyn, ...
                    ];

            fNames = [
                "synodic"
                "adjusted orbital"
                "synodic 2nd"
                ];

        case 'Callisto'

            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    5*fSyn, ...
                    ];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 5th"
                ];

        %% Saturn moons         
        case 'Mimas'
            fTA = 1 / 3600 / 22.67814677274641;
            fOrbTAmid = 1 / 3600 / 22.559683415428385;
            f = [fTA, fOrbTAmid];

            % Harmonics
            f = [f, 2*fOrb, ...
                    ];

            fNames = [
                "true anomaly"
                "adjusted orbital-TA mid"
                "orbital 2nd"
                ];

        case 'Enceladus'
            fTA = 1 / 3600 / 32.927655041360495; % True anomaly period
            f = fTA;

            fNames = "true anomaly";

        case 'Dione'            
            f = fOrbAdj - 2*fSyr;

            fNames = "adjusted orbital-half year beat";

        case 'Rhea'
            f = fOrb;

            fNames = "orbital";

        case 'Titan'
            f = fOrb;

            fNames = "orbital";

        %% Uranus moons   
        case 'Miranda'
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fOrb - fSyn, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    2*fSyn + fOrb, ... % 2nd harmonic beats
                    ];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "orbital-synodic beat"
                "synodic+orbital beat"
                "synodic 2nd+orbital beat"
                ];

        case 'Ariel'
            f = fSyn;

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    ];

            fNames = [
                "synodic"
                "synodic 2nd"
                "synodic 3rd"
                ];

        case 'Umbriel'
            f = fSyn;

            % Harmonics
            f = [f, 2*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    ];

            fNames = [
                "synodic"
                "synodic 2nd"
                "synodic-orbital beat"
                ];

        case 'Titania'
            f = fSyn;

            % Harmonics
            f = [f, 2*fSyn, ...
                    ];

            fNames = [
                "synodic"
                "synodic 2nd"
                ];

        case 'Oberon'
            f = fSyn;

            fNames = "synodic";

        %% Neptune moons
        case 'Triton'
            fSyn = fPar + fOrb; % Retrograde orbit
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, fSyn - fOrb, ... % 1st harmonic beats
                    ];

            fNames = [
                "synodic"
                "orbital"
                "synodic-orbital beat"
                ];
                    
    end
    
    if length(f) ~= length(fNames)
        error(['The number of excitation frequencies does not match the number of names for ' ...
            moonName '. Please update any changed f or fNames lists to match.'])
    end
    [f, indfSort] = sort(f, 'descend');
    fNames = fNames(indfSort);
end
