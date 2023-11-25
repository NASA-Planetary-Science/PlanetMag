function [Brcss, Bpcss, Bzcss] = KS_TailMagShieldCylS3(M, Mode, rmap, pmap, zmap, sphi)
% Get magnetotail shield field in cylindrical System III coordinates.
%
% Parameters
% ----------
% M : int
%   Dimension (MxM) of magnetotail coefficients ``a`` and ``c``, typically ``Mode + 1``.
% Mode : int
%   Parameter selection for magnetotail model. Passing 7 or greater will include all coefficients.
% rmap : double, 1xN
%   Distance from z axis in current sheet coordinates in planetary radii.
% pmap : double, 1xN
%   Azimuthal angle in current sheet coordinates in radians.
% zmap : double, 1xN
%   Distance from current sheet normal plane in planetary radii.
% sphi : double, 1xN
%   East longitude of the plane containing the Sun in System III coordinates in radians.
%
% Returns
% -------
% Brcss, Bpcss, Bzcss : double, 1xN
%   Magnetic field contribution from magnetotail in current sheet cylindrical coordinates in nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Rotate from cylindrical to cartesian in current sheet coordinates
    [xcs, ycs, zcs] = cyl2xyz(rmap, pmap, zmap);
    % MJS note: Khurana's code changes to xJSO, etc. here instead of System III. I think the
    % pseudo- part comes in because we are rotating the xz plane to include the Sun, like in JSO
    % coordinates, but the current sheet cylindrical coordinate system is rotated about the x axis
    % relative to the standard JSO coordinates.
    % Rotate to pJSO coordinates
    [xpJSO, ypJSO, zpJSO] = KS_xyz2pJSO(xcs, ycs, zcs, sphi);
    
    % Get shielded magnetotail field
    [BxpJSO, BypJSO, BzpJSO] = KS_BtailShield(M, Mode, xpJSO, ypJSO, zpJSO);
    % Rotate to System III cartesian
    [Bxcs, Bycs, Bzcs] = KS_BpJSOtoBxyz(BxpJSO, BypJSO, BzpJSO, sphi);
    % Rotate back to cylindrical for output
    [Brcss, Bpcss, Bzcss] = Bxyz2Bcyl(Bxcs, Bycs, Bzcs, pmap);
end
