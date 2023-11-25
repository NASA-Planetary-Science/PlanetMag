function [rho, phi, z] = xyz2cyl(x, y, zin)
% Convert position vectors from cartesian coordinates to cylindrical.
%
% Parameters
% ----------
% x : double, 1xN
%   x coordinate to be converted.
% y : double, 1xN
%   y coordinate to be converted.
% zin : double, 1xN
%   z coordinate to be converted.
%
% Returns
% -------
% rho, phi, z : double, 1xN
%   Cylindrical coordinates of converted position vectors. ``rho`` has the same units as the input
%   coordinates. ``phi`` is in radians.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rho = sqrt(x.^2 + y.^2);
    phi = atan2(y, x);
    z = zin;
end
