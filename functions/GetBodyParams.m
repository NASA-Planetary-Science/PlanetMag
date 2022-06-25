function [RparentEq_km, RparentPol_km, RmoonEq_km, omegaParent, omegaMoon, ...
        Tparent_s, Tmoon_s] = GetBodyParams(moonName)

    if strcmp(moonName,'Io')
        moonID = 501;
        parentID = 599;
    elseif strcmp(moonName,'Europa')
        moonID = 502;
        parentID = 599;
    elseif strcmp(moonName,'Ganymede')
        moonID = 503;
        parentID = 599;
    elseif strcmp(moonName,'Callisto')
        moonID = 504;
        parentID = 599;
    elseif strcmp(moonName,'Mimas')
        moonID = 601;
        parentID = 699;
    elseif strcmp(moonName,'Enceladus')
        moonID = 602;
        parentID = 699;
    elseif strcmp(moonName,'Dione')
        moonID = 604;
        parentID = 699;
    elseif strcmp(moonName,'Rhea')
        moonID = 605;
        parentID = 699;
    elseif strcmp(moonName,'Titan')
        moonID = 606;
        parentID = 699;
    elseif strcmp(moonName,'Miranda')
        moonID = 705;
        parentID = 799;
    elseif strcmp(moonName,'Ariel')
        moonID = 701;
        parentID = 799;
    elseif strcmp(moonName,'Umbriel')
        moonID = 702;
        parentID = 799;
    elseif strcmp(moonName,'Titania')
        moonID = 703;
        parentID = 799;
    elseif strcmp(moonName,'Oberon')
        moonID = 704;
        parentID = 799;
    elseif strcmp(moonName,'Triton')
        moonID = 801;
        parentID = 899;
    else
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
    omegaParent = deg2rad(TparentRotRate(2)) / 86400.0;
    
    % Fetch Europa sidereal rotation rate from PCK
    TmoonRotRate = cspice_bodvcd(moonID, 'PM', 3);
    omegaMoon = deg2rad(TmoonRotRate(2)) / 86400.0;
    
    % Also calculate sidereal periods in seconds
    Tparent_s = 2.0*pi / omegaParent;
    Tmoon_s = 2.0*pi / omegaMoon;

end