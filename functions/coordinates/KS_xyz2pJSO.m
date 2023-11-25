function [xpJSO, ypJSO, zpJSO] = KS_xyz2pJSO(x, y, z, sphi)
% Rotate position vectors from System III cartesian to pseudo-JSO coordinates.
% 
% See KS_BpJSOtoBxyz for a definition of the pseudo-JSO frame.
%
% Parameters
% ----------
% x : double, 1xN
%   x coordinate to be converted in the System III frame.
% y : double, 1xN
%   y coordinate to be converted in the System III frame.
% z : double, 1xN
%   z coordinate to be converted in the System III frame.
% sphi : double, 1xN
%   East longitude of the plane containing the Sun in the System III frame in radians.
%
% Returns
% -------
% xpJSO, ypJSO, zpJSO : double, 1xN
%   Cartesian coordinates of converted vectors in pseudo-JSO frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xpJSO =  x .* cos(sphi) + y .* sin(sphi);
	ypJSO = -x .* sin(sphi) + y .* cos(sphi);
	zpJSO =  z;
end

