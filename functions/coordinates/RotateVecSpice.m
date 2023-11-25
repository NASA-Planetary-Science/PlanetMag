function [VxMoon, VyMoon, VzMoon] = RotateVecSpice(Vx, Vy, Vz, ets, coordSource, coordDest)
% Rotate vectors from one SPICE frame to another.
%
% Parameters
% ----------
% Vx : double, 1xN
%   x-aligned component of vectors to be rotated from the source frame to the destination frame.
% Vy : double, 1xN
%   y-aligned component of vectors to be rotated from the source frame to the destination frame.
% Vz : double, 1xN
%   z-aligned component of vectors to be rotated from the source frame to the destination frame.
% ets : double, 1xN
%   Ephemeris times (ETs) in TDB seconds relative to J2000 for each vector to rotate.
% coordSource : char, 1xC
%   Source frame, with coordinate axes to which the input vectors are referenced.
% coordDest : char, 1xC
%   Destination frame, with coordinate axes to which the output vectors are referenced.
%
% Returns
% -------
% VxMoon, VyMoon, VzMoon : double, 1xN
%   Cartesian components of vectors in the destination frame for each input ET.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Vec = [Vx; Vy; Vz];
    rotMat = cspice_pxform(coordSource, coordDest, ets);
    VecMat(:,1,:) = Vec;
    VecMoon = squeeze(pagemtimes(rotMat, VecMat));
    VxMoon = VecMoon(1,:);
    VyMoon = VecMoon(2,:);
    VzMoon = VecMoon(3,:);
end