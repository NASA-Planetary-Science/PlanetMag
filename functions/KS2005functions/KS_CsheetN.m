function [RNx, RNy, RNz] = KS_CsheetN(xS3, yS3, zS3, localTime_rads, stheta, ctimes, AS_CODED)
% Determine current sheet normal vectors at the given evaluation points for the Khurana and
% Schwarzl (2005) model.
%
% Parameters
% ----------
% xS3 : double, 1xN
%   x component of evaluation point position vector in the System III frame.
% yS3 : double, 1xN
%   y component of evaluation point position vector in the System III frame.
% zS3 : double, 1xN
%   z component of evaluation point position vector in the System III frame.
% localTime_rads : double, 1xN
%   Local time in radians, i.e. east longitude in the JSO frame, which is the azimuthal angle
%   between the plane containing the Sun and the evaluation point.
% stheta : double, 1xN
%   (speculated) Angle between the planet--Sun direction and the evaluation point in radians.
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% RNx, RNy, RNz : double, 1xN
%   Orientation vector of unit normal for the current sheet in System III coordinates at each
%   evaluation point.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delta = 0.1;
    xp = xS3 + delta;
    xm = xS3 - delta;
    yp = yS3 + delta;
    ym = yS3 - delta;
    
    % First calculate x derivatives
    rhop = sqrt(xp.^2 + yS3.^2);
    rhom = sqrt(xm.^2 + yS3.^2);
    phip = atan2(yS3, xp);
    phim = atan2(yS3, xm);
    [xJSOp, yJSOp, ~] = KS_S3CtoJSO(xp, yS3, zS3, ctimes, AS_CODED);
    [xJSOm, yJSOm, ~] = KS_S3CtoJSO(xm, yS3, zS3, ctimes, AS_CODED);
    
    zNS3p = KS_CsheetStruc(rhop, phip, xJSOp, yJSOp, localTime_rads, stheta, AS_CODED);
    zNS3m = KS_CsheetStruc(rhom, phim, xJSOm, yJSOm, localTime_rads, stheta, AS_CODED);
    dzdx = (zNS3p - zNS3m) / (2*delta);
    
    % Now calculate y derivatives
    rhop = sqrt(xS3.^2 + yp.^2);
    rhom = sqrt(xS3.^2 + ym.^2);
    phip = atan2(yp, xS3);
    phim = atan2(ym, xS3);
    [xJSOp, yJSOp, ~] = KS_S3CtoJSO(xS3, yp, zS3, ctimes, AS_CODED);
    [xJSOm, yJSOm, ~] = KS_S3CtoJSO(xS3, ym, zS3, ctimes, AS_CODED);
    
    zNS3p = KS_CsheetStruc(rhop, phip, xJSOp, yJSOp, localTime_rads, stheta, AS_CODED);
    zNS3m = KS_CsheetStruc(rhom, phim, xJSOm, yJSOm, localTime_rads, stheta, AS_CODED);
    dzdy = (zNS3p - zNS3m) / (2*delta);   
        
    RN = sqrt(dzdx.^2 + dzdy.^2 + 1);
    RNx = -dzdx ./ RN;
    RNy = -dzdy ./ RN;
    RNz = 1 ./ RN;
end
