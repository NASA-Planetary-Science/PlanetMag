function [RparentEq_km, RparentPol_km, RmoonEq_km, a_AU, omegaParent_radps, omegaMoon_radps, ...
        Tparent_s, Tmoon_s, nutPrecParent, nutPrecMoon] = GetBodyParams(moonName)

    switch(moonName)
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
        case 'Dione'
            moonID = 604;
            parentID = 699;
        case 'Rhea'
            moonID = 605;
            parentID = 699;
        case 'Titan'
            moonID = 606;
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
    
    nutPrecParent = cspice_bodvcd(floor(parentID/100), 'NUT_PREC_ANGLES', 34);
    % pck00010.tpc does not contain precession info for these bodies:
    nutPrecExcluded = [602, 604];
    if ~any(moonID == nutPrecExcluded)
        nutPrecMoon = cspice_bodvcd(moonID, 'NUT_PREC_PM', 17);
    else
        nutPrecMoon = zeros(8,1);
    end
    
    switch parentID
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