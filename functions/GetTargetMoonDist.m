function r_RM = GetTargetMoonDist(sc, moonName, parentName, ets)
% Get distance between two targets at given ephemeris times.
%
% Uses loaded SPICE kernels to determine the distance from the first target to the second in terms
% of planetary equatorial radii of the second target. Planetary radius is as defined in the
% last-loaded text PCK file. Ephemeris times are passed in terms of seconds relative to J2000.
% Position vectors are evaluated in IAU coordinates for the parent body of the second target
% because these coordinates are always defined in SPICE and we are only finding the line-of-sight
% distance, so the coordinates don't matter.
%
% Parameters
% ----------
% sc : char, 1xC
%   Name of target body (e.g. spacecraft) understood by SPICE when upper-cased.
% moonName : char, 1xD
%   Name of target body from which to determine relative distance. Must match a code name
%   understood by SPICE and present in the loaded kernel pool when upper-cased. Can be any body
%   loaded in the kernel pool for which sufficient data is present to determine sc location
%   relative to this body AND for which a radius can be determined from loaded PCK files.
% parentName : char, 1xE
%   Name of parent body understood by SPICE when upper-cased to use for evaluation of position
%   vectors.
% ets : double, 1xN
%   Ephemeris times in TDB seconds relative to J2000.
%
% Returns
% -------
% r_RM : double, 1xN
%   Distance between target bodies in equatorial radii of the second target body (moonName).

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, ~, RM_km, ~, ~, ~, ~, ~, ~, ~] = GetBodyParams(moonName);
    spkParent = upper(parentName);
    spkS3 = ['IAU_' spkParent];
    [scPos_km, ~] = cspice_spkpos(upper(sc), ets, spkS3, 'NONE', spkParent);
    spkMoon = upper(moonName);
    [spkMoonPos_km, ~] = cspice_spkpos(spkMoon, ets, spkS3, 'NONE', spkParent);
    xDiff = spkMoonPos_km(1,:) - scPos_km(1,:);
    yDiff = spkMoonPos_km(2,:) - scPos_km(2,:);
    zDiff = spkMoonPos_km(3,:) - scPos_km(3,:);
    r_RM = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2) / RM_km;
end