function [r, theta, phi] = xyz2sph(x, y, z)
% Convert position vectors from cartesian coordinates to spherical.
%
% Parameters
% ----------
% x : double, 1xN
%   x coordinate to be converted.
% y : double, 1xN
%   y coordinate to be converted.
% z : double, 1xN
%   z coordinate to be converted.
%
% Returns
% -------
% r, theta, phi : double, 1xN
%   Spherical coordinates of converted position vectors. ``r`` has the same units as the input
%   coordinates. ``theta`` and ``phi`` are in radians.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = sqrt(x.^2 + y.^2 + z.^2);
    theta = acos(z ./ r);
    phi = atan2(y, x);
end
