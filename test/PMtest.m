function PMtest(LIVE_PLOTS, nptsApprox)
% Test script to run each focus feature of P\lanetMag to check for errors.
%
% Parameters
% ----------
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% nptsApprox : int, default=30000
%   Default number of points to use in PlanetMag. This default is much less than that in PlanetMag
%   because this function is primarily intended only for smoke testing.
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
    if ~exist('nptsApprox', 'var'); nptsApprox = 30000; end

    % Clear any loaded SPICE kernels
    cspice_kclear;

    % Remove output files to make sure we generate everything we need within this script
    outDataTxt = dir(fullfile('out/*.txt'));
    nFiles = length(outDataTxt);
    for i=1:nFiles
        delete(fullfile('out', outDataTxt(i).name))
    end
    outDataMat = dir(fullfile('out/*.mat'));
    nFiles = length(outDataMat);
    for i=1:nFiles
        delete(fullfile('out', outDataMat(i).name))
    end
    
    % Earth models
    PlanetMag('Moon', 'Swarm', 'IAU', 1, 1, 1, 1, LIVE_PLOTS, nptsApprox);

    % Jupiter models
    for scName=["Galileo", "Juno"]
        for mName=["Io", "Europa", "Ganymede", "Callisto"]
            moonName = char(mName);
            PlanetMag(moonName, char(scName), 'IAU', 1, 1, 0, 1, LIVE_PLOTS, nptsApprox);
        end
    end
    % Test FFT evaluation/plotting with Europa
    PlanetMag('Europa', 'Galileo', 'IAU', 1, 1, 1, 1, LIVE_PLOTS, nptsApprox);
    PlanetMag('Europa', 'Galileo', 'SPRH', 1, 1, 1, 1, LIVE_PLOTS, nptsApprox);

    % Saturn models
    for mName=["Mimas", "Enceladus", "Dione", "Rhea", "Titan"]
        moonName = char(mName);
        PlanetMag(moonName, 'Voyager', 'IAU', 1, 1, 0, 1, LIVE_PLOTS, nptsApprox);
    end

    % Uranus models
    for mName=["Miranda", "Ariel", "Umbriel", "Titania", "Oberon"]
        moonName = char(mName);
        PlanetMag(moonName, 'Voyager', 'IAU', 1, 1, 0, 1, LIVE_PLOTS, nptsApprox);
    end

    % Neptune models
    PlanetMag('Triton', 'Voyager', 'IAU', 1, 1, 0, 1, LIVE_PLOTS, nptsApprox);

    % Spacecraft data comparison
    CompareEarModels(LIVE_PLOTS);
    CompareJupModels(LIVE_PLOTS);
    CompareSatModels(LIVE_PLOTS);
    CompareUraModels(LIVE_PLOTS);
    CompareNepModels(LIVE_PLOTS);
    
end