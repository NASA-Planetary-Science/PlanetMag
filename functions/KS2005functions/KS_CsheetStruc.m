function zNS3 = KS_CsheetStruc(rho, phi, xJSO, yJSO, localTime_rads, stheta, AS_CODED)
% Get distance of the current sheet from the System III equatorial plane.
%
% Parameters
% ----------
% rho : double, 1xN
%   Distance from z axis of equatorial projection of evaluation point position vector in the 
%   System III frame in planetary radii.
% phi : double, 1xN
%   East longitude of evalaution point in radians.
% xJSO : double, 1xN
%   x-aligned vector component in JSO frame (see KS_S3CtoJSO) in planetary radii.
% yJSO : double, 1xN
%   y-aligned vector component in JSO frame (see KS_S3CtoJSO) in planetary radii.
% localTime_rads : double, 1xN
%   Local time in radians, i.e. east longitude in the JSO frame, which is the azimuthal angle
%   between the plane containing the Sun and the evaluation point.
% stheta : double, 1xN
%   (speculated) Angle between the planet--Sun direction and the evaluation point in radians.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% zNS3 : double, 1xN
%   z distance of current sheet from equatorial plane in the System III frame in planetary radii.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Retrieve current sheet parameterization
    [period, omegaJ_deghr, X0, phip0, obliq, incl_solar, C, psi2, psi4] ...
    = KS_coeffsCsheetStruc(AS_CODED);

    % Avoid zero values in x
    xJSO(abs(xJSO) < 1e-9 & xJSO > 0) =  1e-9;
    xJSO(abs(xJSO) < 1e-9 & xJSO < 0) = -1e-9;
    
    % Calculate transformation parameters
    Alfven = -360 ./ (period * (C(4)*cos(localTime_rads - psi4) + C(5)) );
    delay1 = C(1)*rho + 0.5*C(2)*rho.^2.*cos(localTime_rads - psi2) + 0.5*C(3)*rho.^2;
    delay2 = rho ./ Alfven * omegaJ_deghr * pi/180;
    phip = phip0 - delay1 - delay2;
    poleOffsetMax = obliq + incl_solar;
    rho1 = sqrt((X0*tanh(xJSO/X0)).^2 + yJSO.^2);
    
    % Calculate distance from sheet to equatorial plane
    zNS3 = rho1*tan(poleOffsetMax).*cos(phi - phip) + xJSO.*(1 - tanh(abs(X0./xJSO))).*tan(stheta);
end