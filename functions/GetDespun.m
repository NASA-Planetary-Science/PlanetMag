function [x, y, z] = GetDespun(xyz_Rp, despin)
% Rotate Cartesian coordinates through some angle.
%
% Intended to "de-spin" a rotating coordinate system for convenience in plotting, e.g. in the case
% of a trajectory around a rotating planet that is evaluated in planetocentric coordinates.
%
% Parameters
% ----------
% xyz_Rp : double, 3xN
%   Cartesian position vectors in planetary radii.
% despin : double, 1xN
%   Angles by which to rotate each position vector, prograde (toward :math:`+\hat{\phi}`). Thus, to
%   despin the coordinates by :math:`\Delta\phi`, negate the angles passed.
%
% Returns
% -------
% x, y, z : double, 1xN
%   Rotated Cartesian position vectors in planetary radii.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r_Rp = sqrt(xyz_Rp(1,:).^2 + xyz_Rp(2,:).^2 + xyz_Rp(3,:).^2);
    theta = acos(xyz_Rp(3,:) ./ r_Rp);
    phi = atan2(xyz_Rp(2,:), xyz_Rp(1,:));

    x = r_Rp .* sin(theta) .* cos(phi + despin);
    y = r_Rp .* sin(theta) .* sin(phi + despin);
    z = r_Rp .* cos(theta);
end
