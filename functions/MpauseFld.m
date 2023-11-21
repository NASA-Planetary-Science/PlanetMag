function [mpBvecOut, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, ets, xyz_km, Mdip_nT, Odip_km, ...
    S3coords, parentName, MPmodel, coeffPath, SPHOUT, Nmax)
% Brief description.
%
% Parameters
% ----------
% inputName : type, dims
%   Description.
% xyz_km : double, 3xN
%   Cartesian coordinates of measurement locations in ``S3coords`` frame in km.
% Mdip_nT : double, 1x3
%   Dipole magnetic moment in standard cartesian coordinates, as surface-equivalent nT, i.e. this
%   value times the planet volume yields the magnetic moment in SI units.
% Odip_km : double, 1x3
%   Dipole magnetic moment offset in km from planet barycenter in standard cartesian coordinates.
% MPopt : int
%   Index of magnetopause current magnetic field model. See GetModelOpts for more details and
%   available options.
% SPHOUT : bool, default=0
%   Whether to return vectors aligned to spherical coordinate axes (true) or cartesian (false).
% Nmaxin : int, default=10
%   Maximum degree to which to limit intrinsic field models. Only has an effect if the value passed
%   is less than the lesser of the model maximum degree and the maximum implemented in spherical
%   harmonic calculations (currently 10).
%
% Returns
% -------
% outputName : type, dims
%   Description.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    if ~exist('Nmax', 'var'); Nmax = 10; end

    nHeadLines = 2;
    npts = length(ets);
    
    
    % Get planet-solar-magnetospheric (PSM) coordinates used by some models and in all cases to
    % asses whether we're within the magnetopause
    if strcmp(parentName, 'Saturn')
        parentLetter = 'K';
    else
        parentLetter = parentName(1);
        if strcmp(parentName, 'Uranus')
            % Some Uranus MP models use the Ness et al. (1986) offset, tilted dipole model (OTD)
            % location for the USM origin. See https://doi.org/10.1126/science.233.4759.85.
            Odip_km = 25600 * [-0.02; 0.02; -0.31];
        end
    end
    PSMcoords = [parentLetter 'SM'];
    disp(['Rotating to ' PSMcoords ' coordinates.'])
    rotMatS3_PSM = cspice_pxform(S3coords, PSMcoords, ets);
    xyzMat_km(:,1,:) = xyz_km;
    xyzPSM_km = squeeze(pagemtimes(rotMatS3_PSM, xyzMat_km));
    r_km = sqrt(xyz_km(1,:).^2 + xyz_km(2,:).^2 + xyz_km(3,:).^2);
    
    % Solar wind dynamic pressure equation from AB2005
    pSW_nPa = 2 * 1.67e-27 * (nSW_pcc*1e6) .* (vSW_kms*1e3).^2 * 1e9;
    
    % Dipole magnitude & field at the magnetic equator
    B0_nT = sqrt(Mdip_nT(1)^2 + Mdip_nT(2)^2 + Mdip_nT(3)^2);                
    
    if ~strcmp(MPmodel, 'None')
        
        if strcmp(MPmodel, 'AB2005') || strcmp(MPmodel, 'Bode1994') || contains(MPmodel, ...
                'Engle1992')
            % Distance to magnetopause at subsolar point in planetary radii, Eq. 25 of AB2005
            Rss_Rp = 39.81 ./ pSW_nPa.^(0.23); 
            
            if strcmp(MPmodel, 'AB2005')
                % Alexeev and Belenkaya (2005) magnetopause model for Jupiter:
                % https://doi.org/10.5194/angeo-23-809-2005
                Rp_km = 71400; % as reported in the publication
                xyz_Rp = xyz_km / Rp_km;
                r_Rp = r_km / Rp_km;
                
                MPcoords = [parentLetter 'DSZ']; % Planet-Dipole-Solar-Zenith coordinates
                % Get evaluation position in PDSZ coordinates
                disp(['Rotating to ' MPcoords ' coordinates for ' MPmodel ' magnetopause model.'])
                rotMat = cspice_pxform(S3coords, MPcoords, ets);
                xyzMat_Rp(:,1,:) = xyz_Rp;
                xyzMP_Rp = squeeze(pagemtimes(rotMat, xyzMat_Rp));
                
                rOutDiff_Rp = r_Rp - Rss_Rp;
                if any(rOutDiff_Rp >= 0)
                    warning(['AB2005 magnetopause model is only valid for r < Rss. At least ' ...
                        'one value is greater, with a maximum of ' num2str(max(rOutDiff_Rp)) ...
                        ' Rp greater than Rss.'])
                end
                
                % Get angle of dipole tilt relative to PSM z axis
                MdipPSM_nT = squeeze(pagemtimes(rotMatS3_PSM, Mdip_nT'));
                psi = acos(MdipPSM_nT(3,:) / B0_nT);

                coeffs = dlmread(fullfile(coeffPath, 'coeffsJupiterMpauseAB2005.csv'), ',', ...
                    nHeadLines, 0);
                dPara = coeffs(:,1);
                dPerp = coeffs(:,2);
                Dtail = coeffs(:,3);
                Gtail = coeffs(:,4);

                [mpBr, mpBth, mpBphi] = deal(zeros(1,npts));
                thMP = acos(xyzMP_Rp(3,:) ./ r_Rp);
                phiMP = atan2(xyzMP_Rp(2,:), xyzMP_Rp(1,:));

                disp(['Evaluating ' MPmodel ' magnetopause model.'])
                % Empirical scaling factor resulting from model assumptions in the publication
                mpB0_nT = 1.31 * (100./Rss_Rp).^3 + 3.38 * (100./Rss_Rp).^2;
                for n=1:min([6, Nmax])

                    % Avoid divide-by-zero errors in Bphi
                    thMPadj = thMP;
                    thMPadj(thMPadj == 0) = eps;
                    % Get unnormalized Legendre functions needed for dPara, dPerp
                    Pn0 = LegendreS(n, 0, thMP, 1);
                    Pn1 = LegendreS(n, 1, thMPadj, 1);
                    dPn0 = dLegendreS(n, 0, thMP, 1);
                    dPn1 = dLegendreS(n, 1, thMP, 1);

                    dBr = n*(r_Rp./Rss_Rp).^(n-1) .* ( dPara(n)*sin(psi).*Pn0 ...
                        + dPerp(n)*cos(psi).*cos(phiMP).*Pn1 );
                    dBth =  (r_Rp./Rss_Rp).^(n-1) .* ( dPara(n)*sin(psi).*dPn0 ...
                        + dPerp(n)*cos(psi).*cos(phiMP).*dPn1);
                    dBphi = (r_Rp./Rss_Rp).^(n-1) .* (-dPerp(n)*cos(psi).*sin(phiMP).*Pn1 ...
                        ./ sin(thMPadj));

                    mpBr = mpBr + dBr;
                    mpBth = mpBth + dBth;
                    mpBphi = mpBphi + dBphi;

                end
                mpBr = mpBr .* mpB0_nT;
                mpBth = mpBth .* mpB0_nT;
                mpBphi = mpBphi .* mpB0_nT;
                
            else
                % Engle (1992) magnetopause model for Jupiter: https://doi.org/10.1029/92JA02036
                % with rotation-dependent coefficients as fit by Bode (1994):
                % https://apps.dtic.mil/sti/citations/ADA284857
                Rp_km = 71433.6; % as reported in Bode (1994)
                % Note that a radius is not given in Engle (1992), and as Bode was working with
                % Engle, this value is assumed to have been used in determining the Engle (1992)
                % coeffs.
                xyz_Rp = xyz_km / Rp_km;
                r_Rp = r_km / Rp_km;
                
                % Get dipole moment precession phase angle
                [~, ~, sunLonS3_deg, ~, ~] = GetPosSpice('SUN', parentName, ets/3600);
                dipLon_deg = atan2d(Mdip_nT(2), Mdip_nT(1));
                alpha = deg2rad(mod(sunLonS3_deg - dipLon_deg, 360));
                
                % Fetch model coefficients
                if strcmp(MPmodel, 'Bode1994')
                    % Time-dependent coefficients
                    Gnm = coeffsBode1994G(alpha);
                else
                    % Fixed orientation coefficients
                    GnmStatic = dlmread(fullfile(coeffPath, ['coeffsJupiter' MPmodel 'G.csv']), ...
                        ',', nHeadLines, 0);
                    Gnm = repmat(GnmStatic, 1,1,npts);
                end
                
                % Scale normalization coefficient based from Cn = 13.5 for
                % Rss = 60RJ as assumed in Engle (1992)
                Cn = 13.5 * 60^3 ./ Rss_Rp.^3;
                
                MPcoords = PSMcoords;
                thMP = acos(xyzPSM_km(3,:) / Rp_km ./ r_Rp);
                phiMP = atan2(xyzPSM_km(2,:), xyzPSM_km(1,:));

                disp(['Evaluating ' MPmodel ' magnetopause model.'])
                [mpBr, mpBth, mpBphi] = MPsphericalHarmonic(r_Rp./Rss_Rp, thMP, phiMP, Gnm, ...
                    min([Nmax,10]));
                
                mpBr = mpBr .* Cn;
                mpBth = mpBth .* Cn;
                mpBphi = mpBphi .* Cn;
            end
            
            % Rotate to cartesian for output
            mpBvec = zeros(3,npts);
            [mpBvec(1,:), mpBvec(2,:), mpBvec(3,:)] = Bsph2Bxyz(mpBr, mpBth, mpBphi, thMP, phiMP);
            
            % Use Alexeev and Belenkaya (2005) model to evaluate in/out of
            % magnetopause
            Rss_km = Rss_Rp * Rp_km;
            xMP_km = GetMPsurfAB2005(Rss_km, xyzPSM_km);
            OUTSIDE_MP = xyzPSM_km(1,:) > xMP_km;
            
        elseif strcmp(MPmodel, 'Q3mp') || strcmp(MPmodel, 'SM1996')
            % MP field model described in Schulz and McNab (1996): 
            % https://doi.org/10.1029/95JA02987 . Uses a linear 
            % combination of the coefficients they present for 
            % the psi = 0 and psi = 90 model results.
            if strcmp(MPmodel, 'Q3mp')
                % "Q3mp" magnetopause field model for Uranus as described in
                % Herbert (2009) along with the AH5 model.
                Rss_Rp = 18; % Fixed value used in Herbert (2009) based on Voyager experience
                Rp_km = 25559; % as reported in Herbert (2009)
            elseif strcmp(parentName, 'Neptune')
                % Use values from Masters (2014): https://doi.org/10.1002/2014JA020744
                % for use with generic Schulz and McNab model.
                % Fixed value suggested by Masters (2014) based on synthesizing Voyager results
                Rss_Rp = 25;
                Rp_km = 24765; % as reported in Masters (2014)
            else
                error(['SM1996 behavior for ' parentName 'is not yet defined.'])
            end
            
            Rss_km = Rss_Rp * Rp_km;
            xyz_Rp = xyz_km / Rp_km;
            r_Rp = r_km / Rp_km;
            
            % Get coordinates centered on the dipole z offset
            % in PSM for determining in/out of MP boundary --
            % Strategy for combining psi = 0 and psi = 90 detailed
            % in SM1996 suggests an xy offset is incompatible.
            Odipz_km = [0; 0; Odip_km(3)];
            OdipPSM_km = squeeze(pagemtimes(rotMatS3_PSM, Odipz_km));
            xyzPSMdipO_km = xyzPSM_km - OdipPSM_km;
            % Rotate dipole-centered coordinates to SMP frame, with z axis
            % along the dipole direction, as needed for SM1996 model
            MPcoords = ['SM' parentLetter];
            disp(['Rotating to ' MPcoords ' coordinates.'])
            rotMat = cspice_pxform(PSMcoords, MPcoords, ets);
            xyzPSMdipOmat_km(:,1,:) = xyzPSMdipO_km;
            xyzSMUdipO_Rss = squeeze(pagemtimes(rotMat, xyzPSMdipOmat_km)) / Rss_km;
            r_Rss = r_km / Rss_km;
            thMP = acos(xyzSMUdipO_Rss(3,:) / r_Rss);
            phiMP = atan2(xyzSMUdipO_Rss(2,:), xyzSMUdipO_Rss(1,:));
            
            % Load relative G coefficients
            G_psi0 =  dlmread(fullfile(coeffPath, 'coeffsSM1996psi0G.csv'), ',', nHeadLines, 0);
            G_psi90 = dlmread(fullfile(coeffPath, 'coeffsSM1996psi90G.csv'), ',', nHeadLines, 0);
            % Combine according to the rule described in SM1996 and scale to dipole moment. Negate
            % the expression because SM1996 was Earth-focused, in that they defined the z axis to
            % be opposite the dipole moment.
            Gnm = zeros(10, 11, npts);
            % Get angle of dipole tilt relative to PSM -x axis
            MdipPSM_nT = squeeze(pagemtimes(rotMatS3_PSM, Mdip_nT'));
            psi = cspice_vsep(MdipPSM_nT, repmat([-1;0;0], 1,npts));
            parfor i=1:npts
                Gnm(:,:,i) = -B0_nT/Rss_Rp^3 * (cos(psi(i))*G_psi0 + sin(psi(i))*G_psi90);
            end
            
            disp(['Evaluating ' MPmodel ' magnetopause model.'])
            [mpBr, mpBth, mpBphi] = MPsphericalHarmonic(r_Rss, thMP, phiMP, Gnm, min([Nmax,10]));
            
            % Rotate to cartesian for output
            mpBvec = zeros(3,npts);
            [mpBvec(1,:), mpBvec(2,:), mpBvec(3,:)] = Bsph2Bxyz(mpBr, mpBth, mpBphi, thMP, phiMP);
            
            % Use paraboloid model of Schulz and McNab (1996) as in AH5
            rhoSM1996_km = GetMPsurfSM1996(Rss_km, xyzPSMdipO_km);
            rho_km = sqrt(xyzPSMdipO_km(2,:).^2 + xyzPSMdipO_km(3,:).^2);
            OUTSIDE_MP = rho_km > rhoSM1996_km;
            
        elseif strcmp(MPmodel, 'AE2021')
            % Box harmonic model from Arridge and Eggington (2021) for magnetopause current fields.
            Rss_Rp = 19; % Fixed value inferred by AE2021 from results of Toth et al. (2004)
            Rp_km = 25559;
            Rss_km = Rss_Rp * Rp_km;
            xyz_Rp = xyz_km / Rp_km;
            r_Rp = r_km / Rp_km;
            % Import model coefficients
            aik =  dlmread(fullfile(coeffPath, 'coeffsUranusAE2021a.csv'), ',', nHeadLines, 0);
            bik =  dlmread(fullfile(coeffPath, 'coeffsUranusAE2021b.csv'), ',', nHeadLines, 0);
            cjl =  dlmread(fullfile(coeffPath, 'coeffsUranusAE2021c.csv'), ',', nHeadLines, 0);
            djl =  dlmread(fullfile(coeffPath, 'coeffsUranusAE2021d.csv'), ',', nHeadLines, 0);
            pqrs = dlmread(fullfile(coeffPath, 'coeffsUranusAE2021pqrs.csv'), ',', nHeadLines, 0);
            % Get "attack angle" of dipole vs. upstream direction (approximated by Uranus-Sun
            % direction)
            MdipPSM_nT = squeeze(pagemtimes(rotMatS3_PSM, Mdip_nT'));
            alpha = cspice_vsep(MdipPSM_nT, repmat([1;0;0], 1,npts));
            Psi = pi/2 - alpha;
            
            % AE2021 model uses cartesian USM coordinates
            MPcoords = 'USM';
            xyzPSM_Rp = xyzPSM_km / Rp_km;
            [mpBx, mpBy, mpBz] = MPboxHarmonic(xyzPSM_Rp, Psi, aik, bik, cjl, djl, pqrs);
            mpBvec = [mpBx; mpBy; mpBz];
            
            % Use paraboloid model described in AE2021 based on Shue et al. (1997) model
            xi = 0.72; % Fixed value inferred by AE2021 from Toth et al. (2004)
            thDSZ = acos(xyzPSM_km(1,:) ./ r_km);
            rAE2021_km = GetMPsurfS1997(Rss_km, xi, thDSZ);
            OUTSIDE_MP = r_km > rAE2021_km;
            
        else
            
            error(['MPmodel "' MPmodel '" not recognized.'])
            
        end
    
        disp('Rotating back to S3 coordinates for combining vectors.')
        rotMat = cspice_pxform(MPcoords, S3coords, ets);
        mpBvecMat(:,1,:) = mpBvec;
        mpBvecOut = squeeze(pagemtimes(rotMat, mpBvecMat));
        
        if SPHOUT
            thS3 = acos(xyz_Rp(3,:) ./ r_Rp);
            phiS3 = atan2(xyz_Rp(2,:), xyz_Rp(1,:));
            [mpBvecOut(1,:), mpBvecOut(2,:), mpBvecOut(3,:)] = Bxyz2Bsph(mpBvecOut(1,:), ...
                mpBvecOut(2,:), mpBvecOut(3,:), thS3, phiS3);
        end
        
    else
        % No magnetopause field model; evaluate in/out MP anyways
        
        mpBvecOut = zeros(size(xyz_km));
        % Default subsolar point radii
        switch(parentName)
            case 'Jupiter'
                Rss_Rp = 60;
                Rp_km = 71492;

                % Use Alexeev and Belenkaya (2005) model
                Rss_km = Rss_Rp * Rp_km;
                xMP_km = GetMPsurfAB2005(Rss_km, xyzPSM_km);
                OUTSIDE_MP = xyzPSM_km(1,:) > xMP_km;
                
            case 'Saturn'
                Rp_km = 60268;
                
                % Paraboloid model of Arridge et al. (2006): https://doi.org/10.1029/2005JA011574
                coeffs = dlmread(fullfile(coeffPath, 'coeffsSaturnMpauseAea2006.csv'), ',', ...
                    nHeadLines, 0);
                a = coeffs(1,1:4);
                Rss_km = a(1) * pSW_nPa.^(-a(2)) * Rp_km;
                K = a(3) + a(4) * pSW_nPa;
                thDSZ = acos(xyzPSM_km(1,:) / r_km);
                rAea2006_km = GetMPsurfS1997(Rss_km, K, thDSZ);
                OUTSIDE_MP = r_km > rAea2006_km;
                
            case 'Uranus'
                Rss_Rp = 18;
                Rp_km = 25559;
                
                % Use paraboloid model of Schulz and McNab (1996) as in AH5
                % Get coordinates centered on the dipole z offset                
                xyzDipO_km = xyz_km;
                xyzDipO_km(3,:) = xyzDipO_km(3,:) - Odip_km(3);
                xyzDipOmat_km(:,1,:) = xyzDipO_km;
                xyzPSMdipO_km = squeeze(pagemtimes(rotMatS3_PSM, xyzDipOmat_km));
                Rss_km = Rss_Rp * Rp_km;
                
                rhoSM1996_km = GetMPsurfSM1996(Rss_km, xyzPSMdipO_km);
                rho_km = sqrt(xyzPSMdipO_km(2,:).^2 + xyzPSMdipO_km(3,:).^2);
                OUTSIDE_MP = rho_km > rhoSM1996_km;
                
            case 'Neptune'
                Rss_Rp = 25;
                Rp_km = 24764;
                
                % Use generic Alexeev and Belenkaya (2005) model from
                % Jupiter
                Rss_km = Rss_Rp * Rp_km;
                xMP_km = GetMPsurfAB2005(Rss_km, xyzPSM_km);
                OUTSIDE_MP = xyzPSM_km(1,:) > xMP_km;
                
            otherwise
                error(['parentName "' parentName '" not recognized.'])
        end
        
    end
    
end


function [Br, Bth, Bphi] = MPsphericalHarmonic(r_Rss, thMP, phiMP, Gnm, Nmax)
    % Evaluate external-source spherical harmonics for magnetopause surface
    % currents, i.e. all hnm = 0.
    
    [Br, Bth, Bphi] = deal(zeros(1,length(r_Rss)));
    for n=1:min([10, Nmax])

        % Handle m=0 separately to avoid unneccessary
        % calculations for Bphi, which is always 0 for m=0
        Pn0 =   LegendreS(n, 0, thMP);
        dPn0 = dLegendreS(n, 0, thMP);

        dBr = n*(r_Rss).^(n-1) .* squeeze(Gnm(n,1,:))' .*  Pn0;
        dBth =  (r_Rss).^(n-1) .* squeeze(Gnm(n,1,:))' .* dPn0;
        Br = Br + dBr;
        Bth = Bth + dBth;

        for m=1:n

            % Get Schmidt semi-normalized Legendre functions
            thMPadj = thMP;
            thMPadj(thMPadj == 0) = eps;
            Pnm =   LegendreS(n, m, thMPadj);
            dPnm = dLegendreS(n, m, thMP);

            dBr =    n*(r_Rss).^(n-1) .* squeeze(Gnm(n,m+1,:))' .*  Pnm .* cos(m*phiMP);
            dBth =     (r_Rss).^(n-1) .* squeeze(Gnm(n,m+1,:))' .* dPnm .* cos(m*phiMP);
            dBphi = -m*(r_Rss).^(n-1) .* squeeze(Gnm(n,m+1,:))' .*  Pnm ./ sin(thMPadj) ...
                .* sin(m*phiMP);

            Br = Br + dBr;
            Bth = Bth + dBth;
            Bphi = Bphi + dBphi;

        end
    end
end


function [Bx, By, Bz] = MPboxHarmonic(xyzPSM_Rp, Psi, a, b, c, d, pqrs)
    % Tsyganenko (2002) model (https://doi.org/10.1029/2001JA000219) as formulated by Arridge and
    % Eggington (2021)
    [Bperpx, Bperpy, Bperpz, Bparax, Bparay, Bparaz] ...
        = deal(zeros(1, length(Psi)));
    p = pqrs(1:3, 1);
    q = pqrs(:, 2);
    r = pqrs(1:3, 3);
    s = pqrs(:, 4);
    
    for i=1:3
        for k=1:3
            prScale = sqrt(1/p(i)^2 + 1/r(k)^2);
            abRot = a(i,k)*cos(Psi) + b(i,k)*cos(2*Psi);
            dBperpx = -prScale * abRot .* exp(xyzPSM_Rp(1,:)*prScale) ...
                .* cos(xyzPSM_Rp(2,:)/p(i)) .* sin(xyzPSM_Rp(3,:)/r(i));
            dBperpy =   1/p(i) * abRot .* exp(xyzPSM_Rp(1,:)*prScale) ...
                .* sin(xyzPSM_Rp(2,:)/p(i)) .* sin(xyzPSM_Rp(3,:)/r(i));
            dBperpz =  -1/r(k) * abRot .* exp(xyzPSM_Rp(1,:)*prScale) ...
                .* cos(xyzPSM_Rp(2,:)/p(i)) .* cos(xyzPSM_Rp(3,:)/r(i));
            
            Bperpx = Bperpx + dBperpx;
            Bperpy = Bperpy + dBperpy;
            Bperpz = Bperpz + dBperpz;
        end
    end
    
    for j=1:4
        for l=1:4
            qsScale = sqrt(1/q(j)^2 + 1/s(l)^2);
            cdRot = c(j,l)*sin(Psi) + d(j,l)*sin(2*Psi);
            dBparax = -qsScale * cdRot .* exp(xyzPSM_Rp(1,:)*qsScale) ...
                .* cos(xyzPSM_Rp(2,:)/q(j)) .* cos(xyzPSM_Rp(3,:)/s(l));
            dBparay =   1/q(j) * cdRot .* exp(xyzPSM_Rp(1,:)*qsScale) ...
                .* sin(xyzPSM_Rp(2,:)/q(j)) .* cos(xyzPSM_Rp(3,:)/s(l));
            dBparaz =   1/s(l) * cdRot .* exp(xyzPSM_Rp(1,:)*qsScale) ...
                .* cos(xyzPSM_Rp(2,:)/q(j)) .* sin(xyzPSM_Rp(3,:)/s(l));
            
            Bparax = Bparax + dBparax;
            Bparay = Bparay + dBparay;
            Bparaz = Bparaz + dBparaz;
        end
    end
    
    Bx = Bperpx + Bparax;
    By = Bperpy + Bparay;
    Bz = Bperpz + Bparaz;
end


function rMP_km = GetMPsurfS1997(Rss_km, xi, thDSZ)
    % Shue et al. (1997) model: https://doi.org/10.1029/97JA00196
    rMP_km = Rss_km .* (2 ./ (1 + cos(thDSZ))).^xi;
end


function xMP_km = GetMPsurfAB2005(Rss_km, xyzPSM_km)
    % Alexeev and Belenkaya (2005) model
    xMP_km = Rss_km - (xyzPSM_km(2,:).^2 + xyzPSM_km(3,:).^2)/2 ./ Rss_km;
end


function rhoMP_km = GetMPsurfSM1996(Rss_km, xyzPSMdipO_km)
    % Schulz and McNab (1996) model
    rhoMP_km = 1.7696 .* Rss_km .* ...
        tanh(0.61338./Rss_km .* (Rss_km - xyzPSMdipO_km(1,:))).^(1/2.280);
end
