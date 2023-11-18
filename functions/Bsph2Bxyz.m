function [Bx, By, Bz] = Bsph2Bxyz(Br, Bth, Bphi, theta, phi)
% Convert spherical (r, theta, phi) vector to cartesian (x, y, z) vector.
%
% Convert vector components aligned to spherical coordinates into vector components aligned to
% cartesian axes.
% Source: Arfken, Weber, Harris, Mathematical Methods for Physicists, 7th ed, pg. 199 for the unit
% vectors.
%
% Parameters
% ----------
% Br : double, 1xN
%   Radial (:math:`\hat{r}\upphi\ce{H2SO4^{3-}}`) component of vectors to be converted.
% Bth : double, 1xN
%   Colatitudinal (:math:`\hat{\theta}`) component of vectors to be converted.
% Bphi : double, 1xN
%   Longitudinal (:math:`\hat{\phi}`) component of vectors to be converted.
% theta : double, 1xN
%   Colatitude for each location associated with each (Br, Bth, Bphi) vector in radians.
% phi : double, 1xN
%   East longitude for each location associated with each (Br, Bth, Bphi) vector in radians.
%
% Returns
% -------
% Bx : double, 1xN
%   x-aligned component of converted vectors.
% By : double, 1xN
%   y-aligned component of converted vectors.
% Bz : double, 1xN
%   z-aligned component of converted vectors.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Bx =  sin(theta) .* cos(phi) .* Br ...
        + cos(theta) .* cos(phi) .* Bth ...
        - sin(phi) .* Bphi;
    By =  sin(theta) .* sin(phi) .* Br ...
        + cos(theta) .* sin(phi) .* Bth ...
        + cos(phi) .* Bphi;
    Bz =  cos(theta) .* Br ...
        - sin(theta) .* Bth;
end