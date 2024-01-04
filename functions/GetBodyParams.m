function [RparentEq_km, RparentPol_km, RmoonEq_km, a_AU, omegaParent_radps, omegaMoon_radps, ...
        Tparent_s, Tmoon_s, nutPrecParent, nutPrecMoon] = GetBodyParams(moonName)
% Return physical properties of the target moon, its orbit, and the parent planet.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of the target moon. Currently implemented options are:
%   
%       - ``'Moon'``
%       - ``'Io'``
%       - ``'Europa'``
%       - ``'Ganymede'``
%       - ``'Callisto'``
%       - ``'Mimas'``
%       - ``'Enceladus'``
%       - ``'Dione'``
%       - ``'Rhea'``
%       - ``'Titan'``
%       - ``'Miranda'``
%       - ``'Ariel'``
%       - ``'Umbriel'``
%       - ``'Titania'``
%       - ``'Oberon'``
%       - ``'Triton'``
%
% and their parent planets.
%
% Returns
% -------
% RparentEq_km : double
%   Equatorial radius of the parent planet in km.
% RparentPol_km : double
%   Polar radius of the parent planet in km.
% RmoonEq_km : double
%   Equatorial radius of the target moon in km.
% a_AU : double
%   Semimajor axis of the parent planet's solar orbit in AU.
% omegaParent_radps : double
%   Angular rotation rate of the parent planet in radians/s.
% omegaMoon_radps : double
%   Angular rotation rate of the target moon in radians/s.
% Tparent_s : double
%   Sidereal rotation period of the parent planet in s.
% Tmoon_s : double
%   Sidereal rotation period of the target moon in s.
% nutPrecParent : double, 1xN'
%   Nutation/precession rate coefficients as specified in the latest IAU report implemented in the
%   loaded PCK file. Size depends on the body and which rates are defined.
% nutPrecMoon : double, 1xM'
%   Nutation/precession rate coefficients as specified in the latest IAU report implemented in the
%   loaded PCK file. Size depends on the body and which rates are defined.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch(moonName)
        case 'Moon'
            moonID = 301;
            parentID = 399;
        case 'Io'
            moonID = 501;
            parentID = 599;
        case 'Europa'
            moonID = 502;
            parentID = 599;
        case 'Ganymede'
            moonID = 503;
            parentID = 599;
        case 'Callisto'
            moonID = 504;
            parentID = 599;
        case 'Mimas'
            moonID = 601;
            parentID = 699;
        case 'Enceladus'
            moonID = 602;
            parentID = 699;
        case 'Tethys'
            moonID = 603;
            parentID = 699;
        case 'Dione'
            moonID = 604;
            parentID = 699;
        case 'Rhea'
            moonID = 605;
            parentID = 699;
        case 'Titan'
            moonID = 606;
            parentID = 699;
        case 'Iapetus'
            moonID = 608;
            parentID = 699;
        case 'Miranda'
            moonID = 705;
            parentID = 799;
        case 'Ariel'
            moonID = 701;
            parentID = 799;
        case 'Umbriel'
            moonID = 702;
            parentID = 799;
        case 'Titania'
            moonID = 703;
            parentID = 799;
        case 'Oberon'
            moonID = 704;
            parentID = 799;
        case 'Triton'
            moonID = 801;
            parentID = 899;
        case 'Earth'
            moonID = 399;
            parentID = 399;
        case 'Jupiter'
            moonID = 599;
            parentID = 599;
        case 'Saturn'
            moonID = 699;
            parentID = 699;
        case 'Uranus'
            moonID = 799;
            parentID = 799;
        case 'Neptune'
            moonID = 899;
            parentID = 899;
        otherwise
            error([moonName ' does not have a defined behavior in GetBodyParams.'])
    end
    
    % Fetch parent body radius from PCK file
    RparentTri = cspice_bodvcd(parentID, 'RADII', 3);
    RparentEq_km =  RparentTri(1);
    RparentPol_km = RparentTri(3);
    
    % Fetch triaxial radii for moon
    RmoonTri = cspice_bodvcd(moonID, 'RADII', 3);
    
    % Use mean of polar and equatorial radii for moons
    RmoonEq_km = (RmoonTri(1) + RmoonTri(3)) / 2;

    % Fetch Jupiter period from PCK
    TparentRotRate = cspice_bodvcd(parentID, 'PM', 3);
    omegaParent_radps = deg2rad(TparentRotRate(2)) / 86400.0;
    
    % Fetch Europa sidereal rotation rate from PCK
    TmoonRotRate = cspice_bodvcd(moonID, 'PM', 3);
    omegaMoon_radps = deg2rad(TmoonRotRate(2)) / 86400.0;
    
    % Also calculate sidereal periods in seconds
    Tparent_s = abs(2.0*pi / omegaParent_radps); % abs because Uranus has a negative number here
    Tmoon_s = abs(2.0*pi / omegaMoon_radps); % abs because some moons rotate retrograde
    
    nutPrecParent = zeros(36,1);
    nutPrecParentRet = cspice_bodvcd(floor(parentID/100), 'NUT_PREC_ANGLES', 36);
    nutPrecParent(1:length(nutPrecParentRet)) = nutPrecParentRet;
    % pck00010.tpc does not contain precession info for these bodies:
    nutPrecExcluded = [602, 604];
    nutPrecMoon = zeros(18,1);
    if ~any(moonID == nutPrecExcluded)
        % bodvcd returns a variable length; use max size we might need
        nutPrecMoonRet = cspice_bodvcd(moonID, 'NUT_PREC_PM', 18);
        nutPrecMoon(1:length(nutPrecMoonRet)) = nutPrecMoonRet;
    end
    
    switch parentID
        case 399
            a_AU = 1.000;
        case 599
            a_AU = 5.204;
        case 699
            a_AU = 9.573;
        case 799
            a_AU = 19.165;
        case 899
            a_AU = 30.178;
    end

end
