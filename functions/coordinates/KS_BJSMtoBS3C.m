function [BxS3, ByS3, BzS3] = KS_BJSMtoBS3C(BxJSM, ByJSM, BzJSM, ctimes, AS_CODED)
% Rotate vectors aligned to Jupiter--Sun--magnetic coordinates to System III cartesian.
% 
% In Jupiter--Sun--magnetic (JSM) coordinates, :math:`\hat{x}` points toward the Sun, the
% xz plane contains the magnetic dipole moment, and :math:`\hat{y}` completes the right-handed set.
% This means :math:`\hat{z}` is directed toward the :math:`x=0` plane projection of the magnetic
% dipole moment. For Jupiter, :math:`\hat{y}` nods roughly above and below the direction opposite
% the orbital velocity vector as the planet rotates.
%
% Parameters
% ----------
% BxJSM : double, 1xN
%   x-aligned vector component to be converted in the JSM frame.
% ByJSM : double, 1xN
%   y-aligned vector component to be converted in the JSM frame.
% BzJSM : double, 1xN
%   z-aligned vector component to be converted in the JSM frame.
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% BxS3, ByS3, BzS3 : double, 1xN
%   x-, y-, and z-aligned components of converted vectors in System III frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~AS_CODED
        ets = KS_ctime2et(ctimes);
        % See description of JSM in KS_S3CtoJSM
        rot = cspice_pxform('JSM', 'IAU_JUPITER', ets);
    else
        % Retrieve dipole orientation and initialize rotation matrices
        dipole = KS_coeffsJSMdipole();
        rot = zeros(3,3, length(ctimes));

        % Construct rotation matrix with the JSM basis vectors as rows of the matrix, for the 
        % S3C -> JSM conversion, then transpose to get JSM -> S3C conversion.
        % First, get xJSM basis vector orientation
        [stheta, sphi, ~] = KS_JSun(ctimes, AS_CODED);
        rot(1,1,:) = cos(stheta) .* cos(sphi);
        rot(1,2,:) = cos(stheta) .* sin(sphi);
        rot(1,3,:) = sin(stheta);
        % Next get yJSM as Mdip x xJSM
        rot(2,1,:) = rot(1,3,:) * dipole(3,2) - rot(1,2,:) * dipole(3,3);
        rot(2,2,:) = rot(1,1,:) * dipole(3,3) - rot(1,3,:) * dipole(3,1);
        rot(2,3,:) = rot(1,2,:) * dipole(3,1) - rot(1,1,:) * dipole(3,2);
        denom = sqrt(rot(2,1,:).^2 + rot(2,2,:).^2 + rot(2,3,:).^2);
        rot(2,1,:) = rot(2,1,:) ./ denom;
        rot(2,2,:) = rot(2,2,:) ./ denom;
        rot(2,3,:) = rot(2,3,:) ./ denom;
        % Finally, get zJSM as xJSM x yJSM
        rot(3,1,:) = rot(1,2,:) .* rot(2,3,:) - rot(1,3,:) .* rot(2,2,:);
        rot(3,2,:) = rot(1,3,:) .* rot(2,1,:) - rot(1,1,:) .* rot(2,3,:);
        rot(3,3,:) = rot(1,1,:) .* rot(2,2,:) - rot(1,2,:) .* rot(2,1,:);
        % Now transpose to get S3C -> JSM rotation
        for i=1:length(ctimes)
            rot(:,:,i) = rot(:,:,i)';
        end
    end
    
    BvecIn = [BxJSM; ByJSM; BzJSM];
    BvecOut = zeros(size(BvecIn));
    for i=1:length(ctimes)
        BvecOut(:,i) = rot(:,:,i) * BvecIn(:,i);
    end
    BxS3 = BvecOut(1,:);
    ByS3 = BvecOut(2,:);
    BzS3 = BvecOut(3,:);
end
