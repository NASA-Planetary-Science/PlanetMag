function [Brho, Bphi, Bzout] = Bxyz2Bcyl(Bx, By, Bz, phi)
% Convert cartesian (x, y, z) vector to cylindrical (rho, phi, z) vector.
%
% Convert vector components aligned to cartesian coordinates into vector components aligned to
% cylindrical axes.
% Source: Arfken, Weber, Harris, Mathematical Methods for Physicists, 7th ed, pg. 197 for the unit
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
% phi : double, 1xN
%   East longitude for each location associated with each (Bx, By, Bz) vector in radians.
%
% Returns
% -------
% Brho, Bphi, Bzout : double, 1xN
%   Equatorial projection of radial (:math:`\hat{\rho}`), azimuthal (:math:`\hat{\phi}`), and
%   axial (:math:`\hat{z}`) components of converted vectors.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Brho  =  cos(phi) .* Bx + sin(phi) .* By;
    Bphi  = -sin(phi) .* Bx + cos(phi) .* By;
    Bzout = Bz;
end
