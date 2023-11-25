function [xJSM, yJSM, zJSM] = KS_S3CtoJSM(xS3C, yS3C, zS3c, ctimes, AS_CODED)
% Rotate position vectors from System III cartesian to Jupiter--Sun--magnetic coordinates.
% 
% See KS_BJSMtoBS3C for a definition of the JSM frame.
%
% Parameters
% ----------
% x : double, 1xN
%   x coordinate to be converted in the System III frame.
% y : double, 1xN
%   y coordinate to be converted in the System III frame.
% z : double, 1xN
%   z coordinate to be converted in the System III frame.
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% xJSM, yJSM, zJSM: double, 1xN
%   Cartesian coordinates of converted vectors in JSM frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~AS_CODED
        % MJS note: The juno_v12 frames kernel does not have a definition for the JSM frame. I
        % created frame definitions for the dipole coordinate system (JUNO_JMAG_O4) as described in
        % the following: https://lasp.colorado.edu/home/mop/files/2015/02/CoOrd_systems12.pdf.
        % This document incorrectly describes this frame as equivalent to JSM. I also created a
        % JUNO_JSM frame definition based on the frames kernels required reading and the JSM
        % description at
        % https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=GO-J-POS-6-SC-TRAJ-JUP-COORDS-V1.0
        % The above was later updated to simply 'JSM' in the project custom frame kernel.
        ets = KS_ctime2et(ctimes);
        rot = cspice_pxform('IAU_JUPITER', 'JSM', ets);
    else
        % Retrieve dipole orientation and initialize rotation matrices
        dipole = KS_coeffsJSMdipole();
        rot = zeros(3,3, length(ctimes));

        % Construct rotation matrix with the JSM basis vectors as rows of the matrix
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
    end
    
    posIn = [xS3C; yS3C; zS3c];
    posOut = zeros(size(posIn));
    for i=1:length(ctimes)
        posOut(:,i) = rot(:,:,i) * posIn(:,i);
    end
    xJSM = posOut(1,:);
    yJSM = posOut(2,:);
    zJSM = posOut(3,:);
end
