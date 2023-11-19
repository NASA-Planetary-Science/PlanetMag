function [rMinMoon_km, rPar_km] = GetMinMoonDist(sc, parentName, ets)
% Get distance from a target to the nearest major moon of a planetary system.
%
% Uses loaded SPICE kernels to determine the distance from the target to the major moons of the
% specified parent body at given ephemeris times (ETs). ETs are in terms of seconds relative to
% J2000. Position vectors are evaluated in IAU coordinates for the parent body since we only
% consider line-of-sight distance. Returns the distance to the nearest moon at each ET and the
% radial distance to the parent body center of mass, both in km.
%
% Parameters
% ----------
% sc : char, 1xC
%   Name of target body (e.g. spacecraft) understood by SPICE when upper-cased.
% parentName : char, 1xD
%   Name of parent body from which to determine relative distance to the target and for the moons
%   of this body. Must be one of the following:
%
%       -``Earth``
%       -``Jupiter``
%       -``Saturn``
%       -``Uranus``
%       -``Neptune``
% 
%   and have SPICE kernels present in the loaded kernel pool sufficient to determine all relative
%   locations.
% ets : double, 1xN
%   Ephemeris times in TDB seconds relative to J2000.
%
% Returns
% -------
% rMinMoon_km : double, 1xN
%   Radial distance from the target to the nearest major moon of the parent body in km.
% rPar_km : double, 1xN
%   Radial distance from the target to the parent body in km.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch parentName
        case 'Earth'
            moons = "Moon";
        case 'Jupiter'
            moons = ["Io", "Europa", "Ganymede", "Callisto"];
        case 'Saturn'
            moons = ["Enceladus", "Rhea", "Dione", "Mimas", "Titan"];
        case 'Uranus'
            moons = ["Miranda", "Ariel", "Umbriel", "Titania", "Oberon"];
        case 'Neptune'
            moons = "Triton";
    end
    nMoons = length(moons);
    spkParent = upper(parentName);
    spkS3 = ['IAU_' spkParent];
    [scPos_km, ~] = cspice_spkpos(upper(sc), ets, spkS3, 'NONE', spkParent);
    spkMoonPos_km = zeros(nMoons, 3, length(ets));
    rMoons = zeros(nMoons, length(ets));
    for i=1:nMoons
        spkMoon = upper(char(moons(i)));
        [spkMoonPos_km(i,:,:), ~] = cspice_spkpos(spkMoon, ets, spkS3, 'NONE', spkParent);
        xDiff = squeeze(spkMoonPos_km(i,1,:))' - scPos_km(1,:);
        yDiff = squeeze(spkMoonPos_km(i,2,:))' - scPos_km(2,:);
        zDiff = squeeze(spkMoonPos_km(i,3,:))' - scPos_km(3,:);
        rMoons(i,:) = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2);
    end
    rMinMoon_km = squeeze(min(rMoons, [], 1));
    rPar_km = sqrt(scPos_km(1,:).^2 + scPos_km(2,:).^2 + scPos_km(3,:).^2);
end