function [xJSO, yJSO, zJSO] = KS_S3CtoJSO(xS3C, yS3C, zS3c, ctimes, AS_CODED)
% Rotate position vectors from System III cartesian to Jupiter--Sun--Orbital coordinates.
% 
% In the Jupiter--Sun--Orbital (JSO) frame, :math:`\hat{x}` points toward the Sun, :math:`\hat{y}`
% points in the direction of the **Sun's** velocity vector as observed from the planet center of
% mass (opposite the orbital velocity vector of Jupiter), and :math:`\hat{z}` completes the 
% right-handed set, roughly toward the ecliptic normal.
%
% Parameters
% ----------
% xS3C : double, 1xN
%   x coordinate to be converted in the System III frame.
% yS3C : double, 1xN
%   y coordinate to be converted in the System III frame.
% zS3C : double, 1xN
%   z coordinate to be converted in the System III frame.
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% xJSO, yJSO, zJSO: double, 1xN
%   Cartesian coordinates of converted vectors in JSO frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~AS_CODED
        ets = KS_ctime2et(ctimes);
        rot = cspice_pxform('IAU_JUPITER', 'JSO', ets);
    else
        % Initialize rotation matrices
        rot = zeros(3,3, length(ctimes));
        % First, get xJSO basis vector orientation
        [stheta, sphi, sphase] = KS_JSun(ctimes, AS_CODED);
        rot(1,1,:) = cos(stheta) .* cos(sphi);
        rot(1,2,:) = cos(stheta) .* sin(sphi);
        rot(1,3,:) = sin(stheta);
        % Next, get the zJSO axis from the obliquity
        % Retrieve obliquity angle
        [~, ~, ~, ~, obliq, ~, ~, ~, ~] = KS_coeffsJSun(AS_CODED);
        rot(3,1,:) = sin(obliq) * cos(sphase);
        rot(3,2,:) = sin(obliq) * sin(sphase);
        rot(3,3,:) = cos(obliq);
        % Now get the yJSO axis from Z x X
        rot(2,1,:) = rot(3,2,:) .* rot(1,3,:) - rot(3,3,:) .* rot(1,2,:);
        rot(2,2,:) = rot(3,3,:) .* rot(1,1,:) - rot(3,1,:) .* rot(1,3,:);
        rot(2,3,:) = rot(3,1,:) .* rot(1,2,:) - rot(3,2,:) .* rot(1,1,:);
    end
    
    posIn = [xS3C; yS3C; zS3c];
    posOut = zeros(size(posIn));
    for i=1:length(ctimes)
        posOut(:,i) = rot(:,:,i) * posIn(:,i);
    end
    xJSO = posOut(1,:);
    yJSO = posOut(2,:);
    zJSO = posOut(3,:);
end
