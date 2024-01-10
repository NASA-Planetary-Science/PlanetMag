function PrintSurfaceMap(planet, model, hpDir, hpFname, outDir, outFbase)
% Evaluates magnetic field vectors on a HEALpix map at the body surface and prints to disk.
%
% Prints a .mat file for the given planet and model that contains magnetic field vectors at the
% "surface" of the planet, typically the 1-bar radius defined in IAU reports.
%
% Parameters
% ----------
% planet : char, 1xC
%   Planet name for which to calculate the surface map.
% model : int
%   Model index to pass to GetModelOpts. If -1 is passed, print maps for all models. If 0 is
%   passed, the default model will be evaluated.
% hpDir : char, 1xD, default='publication'
%   Directory containing a file with pixel locations on a HEALpix grid where function outputs are
%   to be compared.
% hpFname : char, 1xE, default='healpix_locs.txt'
%   File name where HEALpix locations are stored.
% outDir : char, 1xF, default='out'
%   Directory where a file with function outputs will be written.
% outFbase : char, 1xG, default='surfaceMap\_'
%   File name pattern for function outputs.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('hpDir', 'var'); hpDir = 'publication'; end
    if ~exist('hpFname', 'var'); hpFname = 'healpix_locs.txt'; end
    if ~exist('outDir', 'var'); outDir = 'publication'; end
    if ~exist('outFbase', 'var'); outFbase = 'surfMap_'; end

    % HEALpix file information
    hpLocsFile = fullfile(hpDir, hpFname);

    % Load HEALpix pixel locations from file
    hpLocs = readtable(hpLocsFile, 'NumHeaderLines', 2);
    theta_rad = hpLocs.Var1';
    phi_rad = hpLocs.Var2';
    varSize = size(theta_rad);
    r_RP = ones(varSize);
    ets = zeros(varSize);
    
    outMatFileBase = fullfile(outDir, outFbase);

    % Retrieve model options
    switch(planet)
        case 'Jupiter'
            opts = 1:7;
        case 'Saturn'
            opts = 1:2;
        case 'Uranus'
            opts = 1:2;
        case 'Neptune'
            opts = 1:1;
    end
    % Overwrite if a specific, valid model number is selected
    if any(ismember(opts, model)) || model == 0
        opts = model:model;
    end

    % Get planet radius
    LoadSpice(planet);
    [Rp_km, ~, ~, ~, ~, ~, ~, ~, ~, ~] = GetBodyParams(planet);
    r_km = r_RP * Rp_km;

    for opt=opts
        [MagModel, CsheetModel, ~, magModelDescrip, fEnd] = GetModelOpts(planet, opt, -1);
        disp(['Evaluating surface map of ' magModelDescrip ' field model for ' planet '.'])
        if contains(magModelDescrip, 'KS2005')
            [Bxyz_nT, ~, ~] = MagFldJupiterKS2005(r_km, theta_rad, phi_rad, ets, 0);
        else
            [Bxyz_nT, ~, ~] = MagFldParent(planet, r_km, theta_rad, phi_rad, ...
                MagModel, CsheetModel, 0, 0);
        end
        Bxyz_G = Bxyz_nT / 1e5;
        Bmag_G = sqrt(Bxyz_G(1,:).^2 + Bxyz_G(2,:).^2 + Bxyz_G(3,:).^2);

        % Rotate to spherical
        [Br, Bth, Bphi] = Bxyz2Bsph(Bxyz_G(1,:), Bxyz_G(2,:), Bxyz_G(3,:), theta_rad, phi_rad);
        Brtp_G = [Br; Bth; Bphi];

        % Save output data
        outFname = [outMatFileBase planet '_' fEnd];
        save(outFname, 'Bxyz_G', 'Brtp_G', 'Bmag_G', 'opt', 'magModelDescrip', '-v7.3');
        disp(['Printed surface map to file: ' outFname '.mat'])
    end

end
