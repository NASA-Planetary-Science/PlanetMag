function [BrS3IMF, BthS3IMF, BphiS3IMF] = KS_GetBIMF(theta, phi, ctimes, AS_CODED)
% Determine IMF field contribution in System III spherical coordinates.
%
% Uses a parameterization from the original Fortran implmenetation of the Khurana and Schwarzl
% (2005) model to evaluate the magnetic field contribution from the interplanetary magnetic field
% (IMF).
%
% Parameters
% ----------
% theta : double, 1xN
%   Colatitude of evaluation points in the System III frame in radians.
% phi : double, 1xN
%   East longitude of evaluation points in the System III frame in radians.
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% BrS3IMF, BthS3IMF, BphiS3IMF : double, 1xN
%   IMF contribution to magnetic field vectors at the evaluation point in System III spherical
%   coordinates in nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [eps1, eps2, ByConst, BzConst] = KS_coeffsBIMF();
    % Get IMF clock angle
    clockAng = atan2(ByConst, BzConst);
    cMultFact = eps1 + eps2*cos(clockAng/2)^2;
    
    % Get IMF in JSM coords
    BxIMF = zeros(size(ctimes));
    ByIMF = cMultFact * ByConst * ones(size(ctimes));
    BzIMF = cMultFact * BzConst * ones(size(ctimes));
    
    % Rotate from JSM to System III cartesian
    [BxS3IMF, ByS3IMF, BzS3IMF] = KS_BJSMtoBS3C(BxIMF, ByIMF, BzIMF, ctimes, AS_CODED);
    % Rotate from cartesian to spherical
    [BrS3IMF, BthS3IMF, BphiS3IMF] = Bxyz2Bsph(BxS3IMF, ByS3IMF, BzS3IMF, theta, phi);
end
