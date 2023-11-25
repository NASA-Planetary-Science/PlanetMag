function [sthetaOut, sphiOut] = KS_GetMappedSunAngle(stheta, sphi, xpx,xpy,xpz, ypx,ypy,ypz, ...
    zpx,zpy,zpz)
% Get angle to Sun in current sheet coordinates.
%
% Parameters
% ----------
% stheta : double, 1xN
%   (speculated) Angle between the planet--Sun direction and the evaluation point in radians.
% sphi : double, 1xN
%   East longitude of the plane containing the Sun in System III coordinates in radians.
% xpx : double, 1xN
%   (speculated) Unit mapping vector x component for current sheet x axis in System III frame.
% xpy : double, 1xN
%   (speculated) Unit mapping vector y component for current sheet x axis in System III frame.
% xpz : double, 1xN
%   (speculated) Unit mapping vector z component for current sheet x axis in System III frame.
% ypx : double, 1xN
%   (speculated) Unit mapping vector x component for current sheet y axis in System III frame.
% ypy : double, 1xN
%   (speculated) Unit mapping vector y component for current sheet y axis in System III frame.
% ypz : double, 1xN
%   (speculated) Unit mapping vector z component for current sheet y axis in System III frame.
% zpx : double, 1xN
%   Current sheet unit normal x component at evaluation point as determined by KS_CsheetN.
% zpy : double, 1xN
%   Current sheet unit normal y component at evaluation point as determined by KS_CsheetN.
% zpz : double, 1xN
%   Current sheet unit normal z component at evaluation point as determined by KS_CsheetN.
%
% Returns
% -------
% sthetaOut : double, 1xN
%   (speculated) Angle between the planet--Sun direction and the evaluation point in radians.
% sphiOut : double, 1xN
%   East longitude of the plane containing the Sun in current sheet coordinates in radians.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [xin, yin, zin] = sph2xyz(1, stheta, sphi);
    xout = xin.*xpx + yin.*xpy + zin.*xpz;
    yout = xin.*ypx + yin.*ypy + zin.*ypz;
    zout = xin.*zpx + yin.*zpy + zin.*zpz;
    [~, sthetaOut, sphiOut] = xyz2sph(xout, yout, zout);
end