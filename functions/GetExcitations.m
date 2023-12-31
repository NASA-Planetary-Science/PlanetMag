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

            f = [f, 2*fSyn];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                ];

        %% Jupiter moons 
        case 'Io'
            fIoTA = 1 / 3600 / 42.315044531808887029; % true anomaly period
            fOrbBeat = 1 / 3600 / 42.289036514276894252;
            fOrbMid1 = 1 / 3600 / 42.431381950208163;
            fOrbMid2 = 1 / 3600 / 42.305626942856122;

            f = [fSyn, fOrb, fIoTA, fOrbMid1, fOrbMid2];

            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    5*fSyn, ...
                    6*fSyn, ...
                    fSyn - fIoTA, ... % 1st harmonic beats
                    fSyn - fOrb, ...
                    fSyn + fOrb, ...
                    fSyn + fIoTA, ...
                    2*fSyn - fIoTA, ... % 2nd harmonic beats
                    2*fSyn + fIoTA, ...
                    3*fSyn - fIoTA, ... % 3rd harmonic beats
                    3*fSyn + fIoTA, ...
                    fSynAdj - fIoTA, ... % 1st harmonic beats
                    fSynAdj - fOrb, ... % 1st harmonic beats
                    fSynAdj + fIoTA, ...
                    2*fSynAdj - fIoTA, ... % 2nd harmonic beats
                    2*fSynAdj + fIoTA, ...
                    3*fSynAdj - fIoTA, ... % 3rd harmonic beats
                    3*fSynAdj + fIoTA, ...
                    fOrbBeat, ... % ~ 2nd harmonic beat between orbital and TA periods
                    fSynAdj - fOrbMid2, ...  % Beat between strong(est) oscillations
                    fPar - fJyr, ... % Solar oscillation periods due to 
                    fPar - 2*fJyr, ... % magnetopause current fields
                    fOrbAdj - fJyr, ...
                    fOrbAdj - 2*fJyr];

            fNames = [
                "synodic"
                "orbital"
                "true anomaly"
                "orbital mid 1"
                "orbital mid 2"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic 5th"
                "synodic 6th"
                "synodic-TA beat"
                "synodic-orbital beat"
                "synodic+orbital beat"
                "synodic+TA beat"
                "synodic 2nd-TA beat"
                "synodic 2nd+TA beat"
                "synodic 3rd-TA beat"
                "synodic 3rd+TA beat"
                "adjusted synodic-TA beat"
                "adjusted synodic-orbital beat"
                "adjusted synodic+TA beat"
                "adjusted synodic 2nd-TA beat"
                "adjusted synodic 2nd+TA beat"
                "adjusted synodic 3rd-TA beat"
                "adjusted synodic 3rd+TA beat"
                "orbital-TA beat 2nd"
                "adjusted synodic-orbital mid 2 beat"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        case 'Europa'
            fTA = 1 / 84.62749 / 3600;
             
            f = [fSyn, fTA, fOrbAdj];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    2*fTA, ...
                    fSyn - fTA, ... % 1st harmonic beats
                    fSyn + fTA, ...
                    2*fSyn - fTA, ... % 2nd harmonic beats
                    2*fSyn + fTA, ...
                    2*fTA - fOrb, ...
                    2*fOrb - fTA, ...
                    fTA + fOrb, ...
                    fTA - fJyr, ... % Solar beats with
                    fOrbAdj - fJyr, ... % orbital periods
                    fPar - fJyr, ... % Solar oscillation periods due to 
                    fPar - 2*fJyr, ... % magnetopause current fields
                    ];

            fNames = [
                "synodic"
                "true anomaly"
                "adjusted orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "TA 2nd"
                "synodic-TA beat"
                "synodic+TA beat"
                "synodic 2nd-TA beat"
                "synodic 2nd+TA beat"
                "TA 2nd-orbital beat"
                "orbital 2nd-TA beat"
                "TA+orbital beat"
                "TA-year beat"
                "adjusted orbital-year beat"
                "planet day-year beat"
                "planet day-half year beat"
                ];

        case 'Ganymede'
            fEuropaTA = 3.282358696657e-6;
            fOrbEuropa = 101.3747235 / 360 / 86400;
            fMystery = 1 / 3600 / 34.724522173543868;
            fMystery2 = 1 / 3600 / 150.271894222476305458;
            f = [fSyn, fOrb];

            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    5*fSyn, ...
                    6*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    2*fSyn - fOrb, ... % 2nd harmonic beats
                    2*fSyn + fOrb, ...
                    3*fSyn - fOrb, ... % 3rd harmonic beats
                    3*fSyn + fOrb, ...
                    fMystery, ...
                    fSyn - fMystery, ...
                    fSyn + fMystery, ...
                    fEuropaTA/2, ... % Half of Europa's true anomaly oscillation
                    fOrbAdj, ...
                    fPar - fOrbEuropa, ...
                    2*fSyn - 5*fOrb, ...
                    fMystery2, ...
                    fSyn - fMystery2, ...
                    1 / 3600 / 4.573237865242181, ...
                    1 / 3600 / 6.208585640926228, ...
                    1 / 3600 / 3.188814749737327, ...
                    1 / 3600 / 3.906250548182006, ...
                    2*fMystery, ...
                    fPar - fJyr, ... % Solar oscillation periods due to 
                    fPar - 2*fJyr, ... % magnetopause current fields
                    fOrbAdj - fJyr, ...
                    fOrbAdj - 2*fJyr];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic 5th"
                "synodic 6th"
                "synodic-orbital beat"
                "synodic+orbital beat"
                "synodic 2nd-orbital beat"
                "synodic 2nd+orbital beat"
                "synodic 3rd-orbital beat"
                "synodic 3rd+orbital beat"
                "mystery"
                "synodic-mystery beat"
                "synodic+mystery beat"
                "Europa TA half period"
                "adjusted orbital"
                "Europa synodic"
                "synodic 2nd-orbital 5th beat"
                "mystery 2"
                "synodic-mystery 2 beat"
                "4.57?"
                "6.21?"
                "3.19?"
                "3.91?"
                "mystery 2nd"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        case 'Callisto'
            fSheet1 = 1 / 3600 / 15.120553206749488;
            fSheet2 = 1 / 3600 / 7.6696338880659605;
            fSheet3 = 1 / 3600 / 3.8072842207111743;
            fSheet4 = 1 / 3600 / 3.059002331397538;

            f = [fSyn, fOrbAdj];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    5*fSyn, ...
                    6*fSyn, ...
                    2*fOrb, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    2*fSyn - fOrb, ... % 2nd harmonic beats
                    2*fSyn + fOrb, ...
                    fSheet1, ...
                    fSheet2, ...
                    fSheet3, ...
                    fSheet4, ...
                    fPar - fJyr, ... % Solar oscillation periods due to 
                    fPar - 2*fJyr, ... % magnetopause current fields
                    fOrbAdj - fJyr, ...
                    fOrbAdj - 2*fJyr];

            fNames = [
                "synodic"
                "adjusted orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic 5th"
                "synodic 6th"
                "orbital 2nd"
                "synodic-orbital beat"
                "synodic+orbital beat"
                "synodic 2nd-orbital beat"
                "synodic 2nd+orbital beat"
                "sheet 1"
                "sheet 2"
                "sheet 3"
                "sheet 4"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        %% Saturn moons         
        case 'Mimas'
            fTA = 1 / 3600 / 22.67814677274641;
            fOrbTAbeat1 = 1 / 3600 / 22.55964638685534;
            fOrbTAbeat2 = 1 / 3600 / 22.50213602844169;
            f = [fOrb, fOrbAdj, fTA, fOrbTAbeat1, fOrbTAbeat2];

            % Harmonics
            f = [f, 2*fOrb, ...  
                    2*fOrbAdj, ...
                    2*fTA, ...  
                    2*fOrbTAbeat1, ...
                    3*fOrb, ...
                    3*fOrbAdj, ...
                    3*fTA, ...
                    fPar - fSyr, ... % Solar oscillation periods due to 
                    fPar - 2*fSyr, ... % magnetopause current fields
                    fOrbAdj - fSyr, ...
                    fOrbAdj - 2*fSyr];

            fNames = [
                "orbital"
                "adjusted orbital"
                "true anomaly"
                "orbital-TA beat 1"
                "orbital-TA beat 2"
                "orbital 2nd"
                "adjusted orbital 2nd"
                "TA 2nd"
                "orbital-TA beat 1 2nd"
                "orbital 3rd"
                "adjusted orbital 3rd"
                "TA 3rd"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        case 'Enceladus'
            fTA = 1 / 3600 / 32.927200612354675; % True anomaly period
            f = [fOrb, fTA];

            f = [f, 2*fTA, ...
                    fPar - fSyr, ... % Solar oscillation periods due to 
                    fPar - 2*fSyr, ... % magnetopause current fields
                    fOrbAdj - fSyr, ...
                    fOrbAdj - 2*fSyr];

            fNames = [
                "orbital"
                "true anomaly"
                "TA 2nd"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        case 'Dione'
            fTA = 1 / 3600 / 65.872263600244693293;
            f = [fOrb, fTA];
            
            f = [f, ...
                    fPar - fSyr, ... % Solar oscillation periods due to 
                    fPar - 2*fSyr, ... % magnetopause current fields
                    fOrbAdj - fSyr, ...
                    fOrbAdj - 2*fSyr];

            fNames = [
                "orbital"
                "true anomaly"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        case 'Rhea'
            f = fOrb;
            
            f = [f, ...
                    fPar - fSyr, ... % Solar oscillation periods due to 
                    fPar - 2*fSyr, ... % magnetopause current fields
                    fOrbAdj - fSyr, ...
                    fOrbAdj - 2*fSyr];

            fNames = [
                "orbital"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        case 'Titan'
            f = fOrb;
            
            f = [f, 2*fOrb, ...
                    fPar - fSyr, ... % Solar oscillation periods due to 
                    fPar - 2*fSyr, ... % magnetopause current fields
                    fOrbAdj - fSyr, ...
                    fOrbAdj - 2*fSyr];

            fNames = [
                "orbital"
                "orbital 2nd"
                "planet day-year beat"
                "planet day-half year beat"
                "adjusted orbital-year beat"
                "adjusted orbital-half year beat"
                ];

        %% Uranus moons   
        case 'Miranda'
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    2*fOrb, ...
                    fOrb - fSyn, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    2*fSyn - fOrb, ... % 2nd harmonic beats
                    2*fSyn + fOrb, ...
                    2*fOrb + fSyn, ...
                    3*fSyn - fOrb, ... % 3rd harmonic beats
                    3*fSyn + fOrb, ...
                    2*fOrb + 2*fSyn]; % Double harmonic beats

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "orbital 2nd"
                "orbital-synodic beat"
                "synodic+orbital beat"
                "synodic 2nd-orbital beat"
                "synodic 2nd+orbital beat"
                "orbital 2nd+synodic beat"
                "synodic 3rd-orbital beat"
                "synodic 3rd+orbital beat"
                "orbital 2nd+synodic 2nd beat"
                ];

        case 'Ariel'
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats % 40.09
                    fSyn + fOrb];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic-orbital beat"
                "synodic+orbital beat"
                ];

        case 'Umbriel'
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic-orbital beat"
                "synodic+orbital beat"
                ];

        case 'Titania'
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic-orbital beat"
                "synodic+orbital beat"
                ];

        case 'Oberon'
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    4*fSyn, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "synodic 4th"
                "synodic-orbital beat"
                "synodic+orbital beat"
                ];

        %% Neptune moons
        case 'Triton'
            fSyn = fPar + fOrb; % Retrograde orbit
            f = [fSyn, fOrb];

            % Harmonics
            f = [f, 2*fSyn, ...
                    3*fSyn, ...
                    2*fOrb, ...
                    fSyn - fOrb, ... % 1st harmonic beats
                    fSyn + fOrb, ...
                    fSyn - 2*fOrb, ...
                    fSyn + 2*fOrb, ...
                    fSyn - 3*fOrb, ...
                    fSyn + 3*fOrb, ...
                    2*fSyn - fOrb, ... % 2nd harmonic beats
                    2*fSyn + fOrb, ...
                    2*fSyn - 2*fOrb, ...
                    2*fSyn + 2*fOrb];

            fNames = [
                "synodic"
                "orbital"
                "synodic 2nd"
                "synodic 3rd"
                "orbital 2nd"
                "synodic-orbital beat"
                "synodic+orbital beat"
                "synodic-orbital 2nd beat"
                "synodic+orbital 2nd beat"
                "synodic-orbital 3rd beat"
                "synodic+orbital 3rd beat"
                "synodic 2nd-orbital beat"
                "synodic 2nd+orbital beat"
                "synodic 2nd-orbital 2nd beat"
                "synodic 2nd+orbital 2nd beat"
                ];
                    
    end
    
    if length(f) ~= length(fNames)
        error(['The number of excitation frequencies does not match the number of names. ' ...
            'Please update any changed f or fNames lists to match.'])
    end
    [f, indfSort] = sort(f, 'descend');
    fNames = fNames(indfSort);
end
