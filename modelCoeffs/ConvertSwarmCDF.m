function ConvertSwarmCDF(fNameIn, fNameOut, fDir)
% Convert magnetic field data files from the Swarm mission from CDF to ASCII text.
%
% Converts CDF file containing Swarm data to .tab text file with positions and field vector
% components aligned to the IAU_EARTH frame. Spice kernels must already be loaded (use
% LoadSpice with ``'Earth', 'Swarm'``). Note that all indices pertaining to timestamps, 
% locations, and measurement data and the corresponding coordinate systems must be set bespoke
% to each CDF file. Converting this function to work with other missions' data products will
% require a careful examination of the second output from cdfread (info) and the ordering of the
% contents in info.Variables.
%
% Parameters
% ----------
% fNameIn : char, 1xC
%   Name of .cdf input file.
% fNameOut : char, 1xD
%   Name of .tab output file.
% fDir : char, 1xE, default=fullfile('MAG', 'Swarm')
%   Directory where input file is located and output file will be printed.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('fDir', 'var'); fDir = fullfile('MAG', 'Swarm'); end
    [data, ~] = cdfread(fullfile(fDir, fNameIn)); % Load CDF file
    [Nmeas, ~] = size(data);

    % Generate SPICE-compatible UTC datetime strings from CDFepoch
    t_UTC = datestr(cellfun(@todatenum, data(:,1)), 'YYYY-mm-DDTHH:MM:SS.FFF');
    
    % Reformat to double arrays
    pos = cell2mat(data(:,3:5));
    latITRF08_deg = pos(:,1);
    lonITRF08_deg = pos(:,2);
    r_km = pos(:,3) / 1e3;

    % Get spherical coordinates
    thITRF08_rad = deg2rad(90 - latITRF08_deg);
    phiITRF08_rad = deg2rad(lonITRF08_deg);

    % Convert from NEC (North-East-Center) coordinates to spherical to use
    % existing transformations
    BNEC_nT = reshape(cell2mat(data(:,12)), [3 Nmeas])';
    Br_nT = -BNEC_nT(:,3); % Center direction is antiradial
    Bth_nT = -BNEC_nT(:,1); % Northward direction is opposite theta-hat direction
    Bphi_nT = BNEC_nT(:,2); % Eastward direction is aligned with phi-hat direction

    % Convert to cartesian so vectors can easily be rotated with spice
    [BxITRF08_nT, ByITRF08_nT, BzITRF08_nT] = Bsph2Bxyz(Br_nT, Bth_nT, Bphi_nT, thITRF08_rad, phiITRF08_rad);

    % Get TDB timestamps
    ets = cspice_str2et(t_UTC);
    % Rotate B from ITRF93 (the only ITRF frame implemented in spice at
    % current, and close to the ITRF08 coordinates used for Swarm) to
    % Earth IAU coordinates
    [Bx_nT, By_nT, Bz_nT] = RotateVecSpice(BxITRF08_nT, ByITRF08_nT, BzITRF08_nT, ets, 'ITRF93', 'IAU_EARTH');
    
    % Now do the same transformation for position vectors
    xITRF08_km = r_km .* sin(thITRF08_rad) .* cos(phiITRF08_rad);
    yITRF08_km = r_km .* sin(thITRF08_rad) .* sin(phiITRF08_rad);
    zITRF08_km = r_km .* cos(thITRF08_rad);
    [xIAU_km, yIAU_km, zIAU_km] = RotateVecSpice(xITRF08_km, yITRF08_km, zITRF08_km, ets, 'ITRF93', 'IAU_EARTH');
    % Get spherical coordinates for conversions and printing
    theta = acos(zIAU_km ./ r_km);
    phi = atan2(xIAU_km, yIAU_km);

    % Convert B to spherical coordinates, which is most convenient for making
    % comparisons in PlanetMag
    [Br_nT, Bth_nT, Bphi_nT] = Bxyz2Bsph(Bx_nT, By_nT, Bz_nT, theta, phi);

    datWidth = 20; % Printed value width including whitespace
    prec = 5; % Precision in digits after decimal for printed values
    Ndat = 6; % Number of data values (field measurements and position)
    
    % Print values to file to match PDS file conventions
    fmtSpec = ['%s' repmat(['%' num2str(datWidth) '.' num2str(5) 'f'], 1, Ndat) '\n'];
    foutID = fopen(fullfile(fDir, fNameOut), 'w');
    for i=1:Nmeas
        fprintf(foutID, fmtSpec, t_UTC(i,:), Br_nT(i), Bth_nT(i), Bphi_nT(i), r_km(i), theta(i), phi(i));
    end
end