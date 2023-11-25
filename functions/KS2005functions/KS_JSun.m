function [stheta, sphi, sphase] = KS_JSun(ctimes, AS_CODED)
% Calculate the direction and phase angle from Jupiter to the Sun as a function of ctime.
%
% Parameters
% ----------
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% stheta : double, 1xN
%   (speculated) Angle between the planet--Sun direction and the evaluation point in radians.
% sphi : double, 1xN
%   East longitude of the plane containing the Sun in System III coordinates in radians.
% sphase : double, 1xN
%   (speculated) Angle between the planet--Sun direction and orbital velocity vector in radians.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize outputs and working variables
    % fphi is phi in Jupiter fixed (non-rotating) coordinates, in degrees
    [stheta, fphi] = deal(zeros(size(ctimes)));

    [yrJup, omegaJ, omegayr, etime1, obliq, tan_ob, ...
        aa, bb, deltaPhi] = KS_coeffsJSun(AS_CODED);
    yrSec = 365.25 * 86400;
    
    t = KS_etimer(ctimes) - etime1;
    x = [cos(omegayr   * t); ...
         sin(omegayr   * t); ...
         cos(2*omegayr * t); ...
         sin(2*omegayr * t); ...
        (t / yrSec).^2; ...
         t / yrSec; ...
         ones(1, length(t))];
     
    for j=1:7
        fphi   =   fphi + bb(j) * x(j,:);
        stheta = stheta + aa(j) * x(j,:);
    end
     
    % Now rotate to Jupiter System III
    % Start with Jupiter orbital motion
    fphi = mod(fphi + t/yrJup*360, 360);
    % Add Jupiter sidereal rotation
    sphi = mod(fphi - t*omegaJ, 360);
    % Convert both angles to radians
    stheta = deg2rad(stheta);
    sphi = deg2rad(sphi);

    % Adjust angles near the magnetotail
    stheta(stheta >= obliq) = obliq;
    % MJS note: I don't think stheta should be able to be negative, so I think the following line
    % is unnecessary.
    stheta(-stheta >= obliq) = -obliq;

    % Now calculate orbital/solar phase (sphase)
    phi21 = mod(sphi + pi + acos(tan(stheta) / tan_ob), 2*pi);
    phi22 = mod(sphi + pi - acos(tan(stheta) / tan_ob), 2*pi);
    dphi = deg2rad(fphi + deltaPhi);
    phi2b = mod(sphi - dphi, 2*pi);

    sphase = phi21;
    dif1 = abs(phi21 - phi2b);
    dif2 = abs(phi22 - phi2b);
    tf = deg2rad(350);
    dif1(dif1 > tf) = 2*pi - dif1(dif1 > tf);
    dif2(dif2 > tf) = 2*pi - dif2(dif2 > tf);
    sphase(dif1 > dif2) = phi22(dif1 > dif2);
end
