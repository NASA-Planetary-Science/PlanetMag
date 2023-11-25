function [Bx, By, Bz] = Bcyl2Bxyz(Brho, Bphi, Bzin, phi)
% Convert cylindrical (rho, phi, z) vector to cartesian (x, y, z).
%
% Convert vector components aligned to cylindrical coordinates into vector components aligned to
% cartesian axes.
% Source: Arfken, Weber, Harris, Mathematical Methods for Physicists, 7th ed, pg. 197 for the unit
% vectors.
%
% Parameters
% ----------
% Brho : double, 1xN
%   x-aligned component of vectors to be converted.
% Bphi : double, 1xN
%   Azimuthal (:math:`\hat{\phi}`) component of vectors to be converted.
% Bzin : double, 1xN
%   Axial (:math:`\hat{z}`) component of vectors to be converted.
% phi : double, 1xN
%   East longitude for each location associated with each (Brho, Bphi, Bzin) vector in radians.
%
% Returns
% -------
% Bx, By, Bz : double, 1xN
%   x-, y-, and z-aligned components of converted vectors.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Bx =  cos(phi) .* Brho - sin(phi) .* Bphi;
    By =  sin(phi) .* Brho + cos(phi) .* Bphi;
    Bz =  Bzin;
end
