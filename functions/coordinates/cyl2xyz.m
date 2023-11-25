function [x, y, z] = cyl2xyz(rho, phi, zin)
% Convert position vectors from cylindrical coordinates to cartesian.
%
% Parameters
% ----------
% rho : double, 1xN
%   rho coordinate to be converted.
% phi : double, 1xN
%   phi coordinate to be converted in radians.
% zin : double, 1xN
%   z coordinate to be converted.
%
% Returns
% -------
% x, y, z : double, 1xN
%   Cartesian coordinates of converted position vectors with the same units as ``rho`` and ``zin``.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = rho .* cos(phi);
    y = rho .* sin(phi);
    z = zin;
end
