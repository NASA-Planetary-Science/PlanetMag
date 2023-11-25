function [x, y, z] = sph2xyz(r, theta, phi)
% Convert position vectors from spherical coordinates to cartesian.
%
% Parameters
% ----------
% r : double, 1xN
%   r coordinate to be converted.
% theta : double, 1xN
%   theta coordinate to be converted in radians.
% phi : double, 1xN
%   phi coordinate to be converted in radians.
%
% Returns
% -------
% x, y, z : double, 1xN
%   Cartesian coordinates of converted position vectors, with units matching ``r``.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
end
