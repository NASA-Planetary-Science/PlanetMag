function [BxpJSO, BypJSO, BzpJSO] = KS_DipoleShielded(parmod, xpJSO, ypJSO, zpJSO, AS_CODED)
% Get shielded dipole field in pseudo-JSO coordinates.
%
% See KS_BpJSOtoBxyz for a definition of the pseudo-Jupiter--Sun--Orbital (pJSO) frame.
%
% Parameters
% ----------
% parmod : struct
%   Must contain fields:
%
%       - **B0** (`double`) -- Magnetic dipole moment in nT in surface-equivalent magnitude.
%       - **dipTilt_deg** (`double`) -- Tilt of dipole moment vector relative to JSM z axis in
%         degrees. This angle varies as Jupiter rotates, but is treated as a constant in the KS2005
%         model.
%       - **Nmodes** (`int`) -- "Number of dipole modes" is how the parameter was labeled in the
%         original Fortran code. Its usage in KS_BMPperp suggests it's some method of indexing
%         Bessel functions and spherical harmonics in that function.
%
% xpJSO : double, 1xN
%   x coordinate of evaluation points in pseudo-JSO frame in planetary radii.
% ypJSO : double, 1xN
%   y coordinate of evaluation points in pseudo-JSO frame in planetary radii.
% zpJSO : double, 1xN
%   z coordinate of evaluation points in pseudo-JSO frame in planetary radii.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% BxpJSO, BypJSO, BzpJSO : double, 1xN
%   Magnetic field contribution from shielded dipole in pseudo-JSO frame aligned with cartesian
%   axes in nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    B0x = parmod.B0 * sind(parmod.dipTilt_deg);
    B0y = 0;
    B0z = parmod.B0 * cosd(parmod.dipTilt_deg);
    % Get dipole field
    [Bxd, Byd, Bzd] = KS_DipoleFld(xpJSO, ypJSO, zpJSO, B0x, B0y, B0z, AS_CODED);
    
    % Get B from magnetopause
    phi = atan2(zpJSO, ypJSO);
    rho = ypJSO .* cos(phi) + zpJSO .* sin(phi);
    [Brho2, Bphi2, Bx2] = KS_BMPperp(rho, phi, xpJSO, parmod.Nmodes);
    
    % Rotate and combine
    By2 = Brho2 .* cos(phi) - Bphi2 .* sin(phi);
    Bz2 = Brho2 .* sin(phi) + Bphi2 .* cos(phi);
    
    BxpJSO = Bxd + Bx2;
    BypJSO = Byd + By2;
    BzpJSO = Bzd + Bz2;
end
