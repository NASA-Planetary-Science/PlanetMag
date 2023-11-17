function [Br, Bth, Bphi] = Bxyz2Bsph(Bx, By, Bz, theta, phi)
% Convert cartesian (x, y, z) vector to spherical (r, theta, phi) vector.
%
% Convert vector components aligned to cartesian axes into vector components aligned to spherical
% coordinates.
% Source: Arfken, Weber, Harris, Mathematical Methods for Physicists, 7th ed, pg. 199 for the unit
% vectors.
%
% Parameters
% ----------
% Bx : double, 1xN
%   x-aligned component of vectors to be converted.
% By : double, 1xN
%   y-aligned component of vectors to be converted.
% Bz : double, 1xN
%   z-aligned component of vectors to be converted.
% theta : double, 1xN
%   Colatitude for each location associated with each (Bx, By, Bz) vector in radians.
% phi : double, 1xN
%   East longitude for each location associated with each (Bx, By, Bz) vector in radians.
%
% Returns
% -------
% Br : double, 1xN
%   Radial (r-hat) component of converted vectors.
% Bth : double, 1xN
%   Colatitudinal (theta-hat) component of converted vectors.
% Bphi : double, 1xN
%   Longitudinal (phi-hat) component of converted vectors.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Br   =  sin(theta) .* cos(phi) .* Bx ...
          + sin(theta) .* sin(phi) .* By ...
          + cos(theta) .* Bz;
    Bth  =  cos(theta) .* cos(phi) .* Bx ...
          + cos(theta) .* sin(phi) .* By ...
          - sin(theta) .* Bz;
    Bphi = -sin(phi) .* Bx ...
          + cos(phi) .* By;
end