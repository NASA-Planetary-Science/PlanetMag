function ptInsideMpause = KS_CheckIfInsideMappedMP(xS3, yS3, zS3, zNS3, ctimes, AS_CODED)
% Return logical array indicating whether each point is inside the mapped magnetopause.
%
% Parameters
% ----------
% xS3 : double, 1xN
%   x coordinate of evaluation points in System III frame in planetary radii.
% yS3 : double, 1xN
%   y coordinate of evaluation points in System III frame in planetary radii.
% zS3 : double, 1xN
%   z coordinate of evaluation points in System III frame in planetary radii.
% zNS3 : double, 1xN
%   z distance of current sheet from equatorial plane in the System III frame in planetary radii.
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% ptInsideMpause : bool, 1xN
%   Whether each evaluation point is inside (true) or outside (false) of the magnetopause.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [xJSM, yJSM, zJSM] = KS_S3CtoJSM(xS3, yS3, zS3+zNS3, ctimes, AS_CODED);
    ptInsideMpause = KS_CheckIfInsideMP(xJSM, yJSM, zJSM);
end
