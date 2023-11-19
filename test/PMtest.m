function PMtest(LIVE_PLOTS)
% Test script to run each focus feature of P\lanetMag to check for errors.
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
%
% Note
% ----
% Setting ``LIVE_PLOTS`` to true will use a lot of system memory, causing a freeze if insufficient
% memory is available.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end

    % Clear any loaded SPICE kernels
    cspice_kclear;

    % Remove output files to make sure we generate everything we need within this script
    outData = dir(fullfile('out/*.txt'));
    nFiles = length(outData);
    for i=1:nFiles
        delete(fullfile('out', outData(i).name))
    end
    
    % Earth models
    PlanetMag('Moon', 'Swarm', 'IAU', 1, 1, 1, 1, LIVE_PLOTS);

    % Jupiter models
    for scName=["Galileo", "Juno"]
        for mName=["Io", "Europa", "Ganymede", "Callisto"]
            moonName = char(mName);
            PlanetMag(moonName, char(scName), 'IAU', 1, 1, 0, 1, LIVE_PLOTS);
        end
    end
    PlanetMag('Europa', 'Galileo', 'IAU', 1, 1, 1, 1, LIVE_PLOTS);
    PlanetMag('Europa', 'Galileo', 'SPRH', 1, 1, 1, 1, LIVE_PLOTS);

    % Saturn models
    for mName=["Mimas", "Enceladus", "Dione", "Rhea", "Titan"]
        moonName = char(mName);
        PlanetMag(moonName, 'Voyager', 'IAU', 1, 1, 0, 1, LIVE_PLOTS);
    end

    % Uranus models
    for mName=["Miranda", "Ariel", "Umbriel", "Titania", "Oberon"]
        moonName = char(mName);
        PlanetMag(moonName, 'Voyager', 'IAU', 1, 1, 0, 1, LIVE_PLOTS);
    end

    % Neptune models
    PlanetMag('Triton', 'Voyager', 'IAU', 1, 1, 0, 1, LIVE_PLOTS);

    % Spacecraft data comparison
    CompareEarModels(LIVE_PLOTS);
    CompareJupModels(LIVE_PLOTS);
    CompareSatModels(LIVE_PLOTS);
    CompareUraModels(LIVE_PLOTS);
    CompareNepModels(LIVE_PLOTS);
    
end