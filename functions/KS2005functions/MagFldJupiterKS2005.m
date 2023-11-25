function [Bvec, Mdip, Odip] = MagFldJupiterKS2005(r_km, theta, phi, ets, SPHOUT, coeffPath, ...
    nHeadLines, AS_CODED)
% Calculate Jupiter magnetic field at evaluation points and times based on the Khurana and Schwarzl
% (2005) model.
%
% Model implementation based on Fortran code from K. Khurana and H. Schwarzl after Khurana and
% Schwarzl (2005) https://doi.org/10.1029/2004JA010757. Much of the original software was
% refactored to create this module, although some functions (e.g. KS_VIP4noDipole, KS_Upot) were
% difficult to interpret and so were translated near-exactly without refactoring. This functions
% are identifiable from their use of non-descript variable naming, i.e. lots of single-letter
% variable names. In several places, more precise model coefficients or corrections to certain
% parameters have been implemented, but can be toggled off by setting ``AS_CODED`` to true, thereby
% matching the behavior of the original Fortran software.
%
% See GetModelOpts or refer to Khurana and Schwarzl (2005) for a description of the model. Note
% that although the embedded magnetopause current model matches the Fortran implementation, most of
% this model was not documented in the code or the publication, so some variable meanings are based
% on speculation.
%
% Parameters
% ----------
% r_km : double, 1xN
%   Radius of evaluation points in km from planet barycenter.
% theta : double, 1xN
%   Planetocentric colatitude of evaluation points in radians. Must be the same length as r_km.
% phi : double, 1xN
%   Planetocentric east longitude of evaluation points in radians. Must be the same length as r_km.
% ets : double, 1xN
%   Ephemeris times associated with evaluation points in TDB seconds relative to J2000. Must be the
%   same length as r_km.
% SPHOUT : bool, default=0
%   Whether to return vectors aligned to spherical coordinate axes (true) or cartesian (false).
% coeffPath : char, 1xE, default='modelCoeffs'
%   Directory where coefficients files are located.
% nHeadLines : int, default=2
%   Number of header lines in model coefficient files. This is typically ``2``: One for a
%   description of the file contents and one for column headers.
% AS_CODED : bool, default=0
%   Whether to use model coefficients and implementation matching the original Fortran code (true)
%   or with increased precision and corrected parameters (false). Specifically, unless ``AS_CODED``
%   is true:
%
%       - Convert distances to units of RJ = 71323 km as in VIP4 paper (Connerney et al., 1998
%         https://doi.org/10.1029/97JA03726)
%       - Use more precise VIP4 tables of g and h
%       - Use corrected values for Jupiter rotation period and rate
%
% Returns
% -------
% Bvec_nT : double, 3xN
%   Magnetic field vector in System III coordinates at evaluation points in nT. Output rows are x,
%   y, z respectively for cartesian or r, theta, phi for spherical.
% Mdip_nT : double, 1x3
%   Dipole magnetic moment in coordinates matching Bvec_nT, as surface-equivalent nT, i.e. this
%   vector times the planet volume yields the magnetic moment in SI units.
% Odip_km : double, 1x3
%   Dipole magnetic moment offset in km from planet barycenter in coordinates matching Bvec_nT.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('nHeadLines', 'var'); nHeadLines = 2; end
    if ~exist('AS_CODED', 'var'); AS_CODED = 0; end
    
    RJ = 71492;    
    r = r_km / RJ;
    ctimes = KS_ctimer(ets);
    npts = length(ets);
    
    % Dimension of magnetotail model coefficients
    M = 8;
    % Mode is the set of parameters to sum over for outer magnetosphere models
    Mode = 7;
    
    % Get VIP4 coefficients in G
    [gVIP4, hVIP4] = KS_coeffsVIP4(coeffPath, nHeadLines, AS_CODED);
    % Convert to nT
    gVIP4 = gVIP4 * 1e5;
    hVIP4 = hVIP4 * 1e5;
    
    % parmod is a list of model parameters.
    if AS_CODED
        parmod.B0 = sqrt(420543.0^2 + 65920.0^2 + 24992.0^2);
    else
        % Magnitude of dipole moment in VIP4 model in nT
        parmod.B0 = sqrt(gVIP4(1,1)^2 + gVIP4(1,2)^2 + hVIP4(1,2)^2);
    end
    % Tilt of dipole axis relative to JSM z axis. This angle oscillates as Jupiter rotates, but the
    % value here is treated as a constant.
    parmod.dipTilt_deg = 0;
    parmod.Nmodes = 3; % "Number of dipole modes" (?)
    
    % Dipole moment vector calculation
    g = gVIP4;
    h = hVIP4;
    M0 = 4*pi*parmod.B0*1e-15*(RJ*1e3)^3 / (4e-7*pi);
    Mdip = [g(1,2), h(1,2), g(1,1)];
    %% Offset dipole -- see Koochak and Fraser-Snith 2017 https://doi.org/10.1002/2017EA000280
    % Note two typos in Koochak and Fraser-Snith (2017) Eq. 4: G11 and G20 should be g11 and g20,
    % and the g20 inside square brackets should be g21 -- see Fraser-Snith (1987)
    % https://doi.org/10.1029/RG025i001p00001.
    % L0 = 2*g10*g20 + sqrt(3)*(g11*g21 + h11*h21)
    % L1 =  -g11*g20 + sqrt(3)*(g10*g21 + g11*g22 + h11*h22)
    % L2 =  -h11*g20 + sqrt(3)*(g10*h21 - h11*g22 + g11*h22)
    % E = (L0*g10 + L1*g11 + L2*h11) / (4*B0^2)
    L0 = 2*g(1,1)*g(2,1) + sqrt(3)*(g(2,1)*g(2,2) + h(1,2)*h(2,2));
    L1 =  -g(1,2)*g(2,1) + sqrt(3)*(g(1,1)*g(2,2) + g(1,2)*g(2,3) + h(1,2)*h(2,3));
    L2 =  -h(1,2)*g(2,1) + sqrt(3)*(g(1,1)*h(2,2) - h(1,2)*g(2,3) + g(1,2)*h(2,3));
    E = (L0*g(1,1) + L1*g(1,2) + L2*h(1,2)) / (4*parmod.B0^2);

    % Unitless offset parameters
    % x-axis offset: eta =  (L1 - g11*E) / (3*B0^2)
    % y-axis offset: zeta = (L2 - h11*E) / (3*B0^2)
    % z-axis offset: xi =   (L0 - g10*E) / (3*B0^2)
    xi =   (L0 - g(1,1)*E) / (3*parmod.B0^2);
    eta =  (L1 - g(1,2)*E) / (3*parmod.B0^2);
    zeta = (L2 - h(1,2)*E) / (3*parmod.B0^2);

    % Dipole offset in km
    if SPHOUT
        rO_km = sqrt(eta^2 + zeta^2 + xi^2);
        thO_rad = acos(xi / rO_km);
        phiO_rad = atan2(zeta, eta);
        Odip = RJ * [rO_km, thO_rad, phiO_rad];
    else
        Odip = RJ * [eta, zeta, xi];
    end
    
    rho = r .* sin(theta);
    xS3 = rho .* cos(phi);
    yS3 = rho .* sin(phi);
    zS3 = r .* cos(theta);
    
    % Convert Jupiter System III cartesian to Jupiter--Sun--Orbital coordinates
    [xJSO, yJSO, zJSO] = KS_S3CtoJSO(xS3, yS3, zS3, ctimes, AS_CODED);
    % Get local solar time
    localTime_rads = atan2(yJSO, xJSO);
    
    % Get Sun angles
    [stheta, sphi, sphase] = KS_JSun(ctimes, AS_CODED);

    % Get unit normals of zp axis
    [RNx, RNy, RNz] = KS_CsheetN(xS3, yS3, zS3, localTime_rads, stheta, ctimes, AS_CODED);
    
    %% Calculate mapped locations in current sheet coordinates
    % z axis:
    zpx = RNx;
    zpy = RNy;
    zpz = RNz;
    zp = xS3.*zpx + yS3.*zpy + zS3.*zpz;
    % y axis:
    yp1 = sqrt(zpx.^2 + zpy.^2);
    ypx =  zpy ./ yp1;
    ypy = -zpx ./ yp1;
    ypz = zeros(1,npts);
    yp = xS3.*ypx + yS3.*ypy + zS3.*ypz;
    % x axis:
    xpx = ypy.*zpz - ypz.*zpy;
    xpy = ypz.*zpx - ypx.*zpz;
    xpz = ypx.*zpy - ypy.*zpx;
    xp = xS3.*xpx + yS3.*xpy + zS3.*xpz;
    
    % Get derivative grid
	dxpdx = xpx;
	dypdx = ypx;
	dzpdx = zpx;
	dxpdy = xpy;
	dypdy = ypy;
	dzpdy = zpy;
	dxpdz = xpz;
	dypdz = ypz;
	dzpdz = zpz;
        
    % Calculate T matrix
    Txx = (dypdy .* dzpdz - dypdz .* dzpdy);
	Txy = (dxpdz .* dzpdy - dxpdy .* dzpdz);
	Txz = (dxpdy .* dypdz - dxpdz .* dypdy);
	Tyx = (dypdz .* dzpdx - dypdx .* dzpdz);
	Tyy = (dxpdx .* dzpdz - dxpdz .* dzpdx);
	Tyz = (dxpdz .* dypdx - dxpdx .* dypdz);
	Tzx = (dypdx .* dzpdy - dypdy .* dzpdx);
	Tzy = (dxpdy .* dzpdx - dxpdx .* dzpdy);
	Tzz = (dxpdx .* dypdy - dxpdy .* dypdx);
    
    % Get distance from System III equatorial plane to current sheet
    zNS3 = KS_CsheetStruc(rho, phi, xJSO, yJSO, localTime_rads, stheta, AS_CODED);
    % Get mapped coordinates in current sheet frame
    [rmap, pmap, zmap] = xyz2cyl(xp, yp, zS3 - zNS3);
    
    %% Calculate field at mapped location in cylindrical coordinates
    scol = pi/2 - stheta;
    [scolOut, sphiOut] = KS_GetMappedSunAngle(scol, sphi, xpx,xpy,xpz, ypx,ypy,ypz, zpx,zpy,zpz);
        
    [Brds, Bpds, Bzds] = KS_DipoleShieldCylS3(parmod, rmap, pmap, zmap, sphiOut, AS_CODED);
    [Brcs, Bpcs, Bzcs] = KS_TailMagNoTilt(Mode, rmap, pmap, zmap, localTime_rads);
    [Brcss, Bpcss, Bzcss] = KS_TailMagShieldCylS3(M, Mode, rmap, pmap, zmap, sphiOut);
        
    % Sum each contribution
    Brmap = Brds + Brcs + Brcss;
    Bpmap = Bpds + Bpcs + Bpcss;
    Bzmap = Bzds + Bzcs + Bzcss;
    
    %% Rotate to cartesian to apply T matrix
    [BXcarMap, BYcarMap, BZcarMap] = Bcyl2Bxyz(Brmap, Bpmap, Bzmap, pmap);
    
    Bxfinal = Txx.*BXcarMap + Txy.*BYcarMap + Txz.*BZcarMap;
    Byfinal = Tyx.*BXcarMap + Tyy.*BYcarMap + Tyz.*BZcarMap;
    Bzfinal = Tzx.*BXcarMap + Tzy.*BYcarMap + Tzz.*BZcarMap;
    
    %% Rotate to spherical coordinates to sum VIP4 model vectors
    [Brfinal, Bthfinal, Bphifinal] = Bxyz2Bsph(Bxfinal, Byfinal, Bzfinal, theta, phi);
    
    [BrVIP4, BthVIP4, BphiVIP4] = KS_VIP4noDipole(r, theta, phi, gVIP4, hVIP4, AS_CODED);
    Br = Brfinal + BrVIP4;
    Bth = Bthfinal + BthVIP4;
    Bphi = Bphifinal + BphiVIP4;
    
    %% Adjust for whether each point is inside magnetopause or not
    [BrS3IMF, BthS3IMF, BphiS3IMF] = KS_GetBIMF(theta, phi, ctimes, AS_CODED);
    ptInsideMpause = KS_CheckIfInsideMappedMP(xS3, yS3, zS3, zNS3, ctimes, AS_CODED);
    
    Br(ptInsideMpause) = Br(ptInsideMpause) + BrS3IMF(ptInsideMpause);
    Bth(ptInsideMpause) = Bth(ptInsideMpause) + BthS3IMF(ptInsideMpause);
    Bphi(ptInsideMpause) = Bphi(ptInsideMpause) + BphiS3IMF(ptInsideMpause);
    Br(~ptInsideMpause) = BrS3IMF(~ptInsideMpause);
    Bth(~ptInsideMpause) = BthS3IMF(~ptInsideMpause);
    Bphi(~ptInsideMpause) = BphiS3IMF(~ptInsideMpause);
    
    % Convert field vectors to cartesian for output
    if SPHOUT
        Bvec = [Br; Bth; Bphi];
    else
        [Bx, By, Bz] = Bsph2Bxyz(Br, Bth, Bphi, theta, phi);
        Bvec = [Bx; By; Bz];
    end
end
