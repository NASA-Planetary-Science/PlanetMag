function [r_km, theta_rad, phi_rad, xyz_km, S3coords] = GetPosSpice(moonName, parentName, t_h, ...
    S3coords)
% Retrieve locations of a target relative to a parent body at specified ephemeris times.
%
% Uses loaded SPICE kernels to determine the position of one body relative to another in a specific
% coordinate system. If the coordinates are not passed, standard cartesian coordinates are
% selected. Ephemeris times are passed in terms of hours relative to J2000.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of target body, which must match a code name understood by SPICE and present in the loaded
%   kernel pool when upper-cased. Can be any body loaded in the kernel pool for which sufficient
%   data is present to determine location relative to the parent.
% parentName : char, 1xD
%   Name of parent body understood by SPICE when upper-cased.
% t_h : double, 1xN
%   Ephemeris times in TDB hours relative to J2000.
% S3coords : char, 1xE, optional
%   Standard (System III) coordinates to use for the parent body. Defaults to ULS for Uranus, NLS
%   for Neptune, or IAU_PARENT for any other parent body. Any coordinate system understood by SPICE
%   will return a result, but the first 3 return values from this function will only be converted
%   correctly into spherical coordinates if cartesian coordinates are specified. Magnetic field
%   model evaluations will only work correctly with the default coordinates.
%
% Returns
% -------
% r_km : double, 1xN
%   Radii from parent body center of mass in km.
% theta_rad : double, 1xN
%   Colatitude of positions relative to parent body spin pole in radians.
% phi_rad : double, 1xN
%   East longitude of positions relative to IAU-defined parent body prime meridian in radians.
% xyz_km : double, 3xN
%   Position vectors of target body for each time in the selected coordinates.
% S3coords : char, 1xE
%   Coordinate system of position vectors.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = t_h * 3600;
    spkParent = upper(parentName);
    spkMoon = upper(moonName);
    if ~exist('S3coords', 'var')
        switch(parentName)
            case 'Uranus'
                S3coords = 'ULS';
            case 'Neptune'
                S3coords = 'NLS';
            otherwise
                S3coords = ['IAU_' spkParent];
        end
    end
    [xyz_km, ~] = cspice_spkpos(spkMoon, t, S3coords, 'NONE', spkParent);
    r_km = sqrt(xyz_km(1,:).^2 + xyz_km(2,:).^2 + xyz_km(3,:).^2);
    theta_rad = acos(xyz_km(3,:) ./ r_km);
    phi_rad = atan2(xyz_km(2,:), xyz_km(1,:));
end