function [Brds, Bpds, Bzds] = KS_DipoleShieldCylS3(parmod, rmap, pmap, zmapin, sphi, AS_CODED)
% Get dipole shield field in cylindrical System III coordinates.
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
% rmap : double, 1xN
%   Distance from z axis in current sheet coordinates in planetary radii.
% pmap : double, 1xN
%   Azimuthal angle in current sheet coordinates in radians.
% zmapin : double, 1xN
%   Distance from current sheet normal plane in planetary radii.
% sphi : double, 1xN
%   East longitude of the plane containing the Sun in System III coordinates in radians.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% Brds, Bpds, Bzds : double, 1xN
%   Magnetic field contribution from shielded dipole aligned to cylindrical System III axes in nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get mapped cartesian coords
    [xmap, ymap, zmap] = cyl2xyz(rmap, pmap, zmapin);
    % Transform to pseudo-JSO position
    [xpJSO, ypJSO, zpJSO] = KS_xyz2pJSO(xmap, ymap, zmap, sphi);
    % Get shielded dipole in pseudo-JSO coordinates
    [BxpJSO, BypJSO, BzpJSO] = KS_DipoleShielded(parmod, xpJSO, ypJSO, zpJSO, AS_CODED);
    % Transform psuedo-JSO field vectors to JS3
    [Bx, By, Bz] = KS_BpJSOtoBxyz(BxpJSO, BypJSO, BzpJSO, sphi);
    % Transform to cylindrical for output
    [Brds, Bpds, Bzds] = Bxyz2Bcyl(Bx, By, Bz, pmap);
end
