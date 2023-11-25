function [Bxd, Byd, Bzd] = KS_DipoleFld(x, y, z, B0x, B0y, B0z, AS_CODED)
% Calculate dipole field as a function of position and dipole components
%
% Parameters
% ----------
% x : double, 1xN
%   x coordinate of evalaution points in the System III frame.
% y : double, 1xN
%   y coordinate of evalaution points in the System III frame.
% z : double, 1xN
%   z coordinate of evalaution points in the System III frame.
% B0x : double, 1xN
%   x-aligned component of dipole moment vector in surface-equivalent magnitude in the System III
%   frame in nT.
% B0y : double, 1xN
%   y-aligned component of dipole moment vector in surface-equivalent magnitude in the System III
%   frame in nT.
% B0z : double, 1xN
%   z-aligned component of dipole moment vector in surface-equivalent magnitude in the System III
%   frame in nT.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% Bxd, Byd, Bzd : double, 1xN
%   Magnetic field contribution from magnetic dipole moment aligned to cartesian System III axes in
%   nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~AS_CODED
        % Convert distances to RJ_VIP4 to fit with that model
        conv = 71492/71323; x = x * conv; y = y * conv; z = z * conv;
    end
    
    r = sqrt(x.^2 + y.^2 + z.^2);
	a11 = (3*x.^2 - r.^2) ./ r.^5;
	a12 = (3*x.*y) ./ r.^5;
	a13 = (3*x.*z) ./ r.^5;
	a21 = (3*x.*y) ./ r.^5;
	a22 = (3*y.^2 - r.^2) ./ r.^5;
	a23 = (3*y.*z) ./ r.^5;
	a31 = (3*x.*z) ./ r.^5;
	a32 = (3*y.*z) ./ r.^5;
	a33 = (3*z.^2 - r.^2) ./ r.^5;
	Bxd = a11.*B0x + a12.*B0y + a13.*B0z;
	Byd = a21.*B0x + a22.*B0y + a23.*B0z;
	Bzd = a31.*B0x + a32.*B0y + a33.*B0z;
end
