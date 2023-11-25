function [Bx, By, Bz] = KS_BpJSOtoBxyz(BxpJSO, BypJSO, BzpJSO, sphi)
% Rotate cartesian vector aligned with pseudo-JSO coordinates to System III.
% 
% In psuedo-Jupiter--Sun--Orbital (pJSO) coordinates, :math:`\hat{z}` is parallel to the jovian
% spin axis and so :math:`\hat{x}` does not point toward the Sun as in JSO, but the Sun does lie in
% the xz plane. pJSO is equivalent to the Jupiter--De-Spun--Sun frame (JUNO_JSS) in the Juno frames
% kernel. See KS_S3CtoJSO for a definition of the JSO frame.
%
% Parameters
% ----------
% BxpJSO : double, 1xN
%   x-aligned vector component to be converted in the pseudo-JSO frame.
% BypJSO : double, 1xN
%   y-aligned vector component to be converted in the pseudo-JSO frame.
% BzpJSO : double, 1xN
%   z-aligned vector component to be converted in the pseudo-JSO frame.
% sphi : double, 1xN
%   East longitude of the plane containing the Sun in the System III frame in radians.
%
% Returns
% -------
% Bx, By, Bz : double, 1xN
%   x-, y-, and z-aligned components of converted vectors in System III frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Bx = BxpJSO .* cos(sphi) - BypJSO .* sin(sphi);
    By = BxpJSO .* sin(sphi) + BypJSO .* cos(sphi);
    Bz = BzpJSO;
end
