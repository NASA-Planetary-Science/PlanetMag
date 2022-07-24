function [mpBvecOut, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, ets, xyz_km, Mdip_nT, parentName, MPmodel, SPHOUT, Nmax)
    coeffPath = './modelCoeffs/';
    nHeadLines = 2;
    npts = length(ets);
    S3coords = ['IAU_' upper(parentName)];
    
    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    if ~exist('Nmax', 'var'); Nmax = 10; end
    
    % Get planet-solar-magnetospheric (PSM) coordinates used by
    % some models and in all cases to asses whether we're within
    % the magnetopause
    if strcmp(parentName, 'Saturn')
        parentLetter = 'K';
    else
        parentLetter = parentName(1);
    end
    PSMcoords = [parentLetter 'SM'];
    xyzPSM_km = zeros(3,npts);
    disp(['Rotating to ' PSMcoords ' coordinates.'])
    rotMat = cspice_pxform(S3coords, PSMcoords, ets);
    parfor i=1:npts
        xyzPSM_km(:,i) = rotMat(:,:,i) * xyz_km(:,i);
    end
    
    % Solar wind dynamic pressure equation from AB2005
    pSW_nPa = 2 * 1.67e-27 * (nSW_pcc*1e6) .* (vSW_kms*1e3).^2 * 1e9;
    
    if ~strcmp(MPmodel, 'None')
        
        if strcmp(MPmodel, 'AB2005') || strcmp(MPmodel, 'Bode1994') || strcmp(MPmodel(1:9), 'Engle1992')
            % Distance to magnetopause at subsolar point in planetary
            % radii, Eq. 25 of AB2005
            Rss_Rp = 39.81 ./ pSW_nPa.^(0.23); 
            
            if strcmp(MPmodel, 'AB2005')
                % Alexeev and Belenkaya (2005) magnetopause model for Jupiter:
                % https://doi.org/10.5194/angeo-23-809-2005
                Rp_km = 71400; % as reported in the publication
                xyz_Rp = xyz_km / Rp_km;
                r_Rp = sqrt(xyz_Rp(1,:).^2 + xyz_Rp(2,:).^2 + xyz_Rp(3,:).^2);

                MPcoords = [parentLetter 'DSZ']; % Planet-Dipole-Solar-Zenith coordinates
                % Get evaluation position in PDSZ coordinates
                xyzMP_Rp = zeros(3,npts);
                disp(['Rotating to ' MPcoords ' coordinates for ' MPmodel ' magnetopause model.'])
                rotMat = cspice_pxform(S3coords, MPcoords, ets);

                parfor i=1:npts
                    xyzMP_Rp(:,i) = rotMat(:,:,i) * xyz_Rp(:,i);
                end
                B0_nT = sqrt(Mdip_nT(1)^2 + Mdip_nT(2)^2 + Mdip_nT(3)^2); % Dipole magnitude & field at the magnetic equator
                rOutDiff_Rp = r_Rp - Rss_Rp;
                if any(rOutDiff_Rp >= 0)
                    warning(['AB2005 magnetopause model is only valid for r < Rss. At least one value is' ...
                        'greater, with a maximum of ' num2str(max(rOutDiff_Rp)) ' Rp greater than Rss.'])
                end

                psi = acos(Mdip_nT(3) / B0_nT); % Angle of dipole tilt relative to spin axis

                coeffs = dlmread(fullfile([coeffPath 'coeffsJupiterMpauseAB2005.csv']), ',', nHeadLines, 0);
                dPara = coeffs(:,1);
                dPerp = coeffs(:,2);
                Dtail = coeffs(:,3);
                Gtail = coeffs(:,4);

                [mpBr, mpBth, mpBphi] = deal(zeros(1,npts));
                thMP = acos(xyzMP_Rp(3,:) ./ r_Rp);
                phiMP = atan2(xyzMP_Rp(2,:), xyzMP_Rp(1,:));

                disp(['Evaluating ' MPmodel ' magnetopause model.'])
                % Empirical scaling factor resulting from model assumptions
                % in the publication
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

                    dBr = n*(r_Rp./Rss_Rp).^(n-1) .* ( dPara(n)*sin(psi)*Pn0  + dPerp(n)*cos(psi)*cos(phiMP).*Pn1 );
                    dBth =  (r_Rp./Rss_Rp).^(n-1) .* ( dPara(n)*sin(psi)*dPn0 + dPerp(n)*cos(psi)*cos(phiMP).*dPn1);
                    dBphi = (r_Rp./Rss_Rp).^(n-1) .* (-dPerp(n)*cos(psi)*sin(phiMP).*Pn1./sin(thMPadj));

                    mpBr = mpBr + dBr;
                    mpBth = mpBth + dBth;
                    mpBphi = mpBphi + dBphi;

                end
                mpBr = mpBr .* mpB0_nT;
                mpBth = mpBth .* mpB0_nT;
                mpBphi = mpBphi .* mpB0_nT;
                
            else
                % Engle (1992) magnetopause model for Jupiter:
                % https://doi.org/10.1029/92JA02036
                % with rotation-dependent coefficients as fit by Bode
                % (1994): https://apps.dtic.mil/sti/citations/ADA284857
                Rp_km = 71433.6; % as reported in Bode (1994)
                % Note that a radius is not given in Engle (1992), and as
                % Bode was working with Engle, this value is assumed to
                % have been used in determining the Engle (1992) coeffs.
                xyz_Rp = xyz_km / Rp_km;
                r_Rp = sqrt(xyz_Rp(1,:).^2 + xyz_Rp(2,:).^2 + xyz_Rp(3,:).^2);
                
                % Get dipole moment precession phase angle
                [~, ~, sunLonS3_deg, ~] = GetPosSpice('SUN', parentName, ets/3600);
                dipLon_deg = atan2d(Mdip_nT(2), Mdip_nT(1));
                alpha = deg2rad(mod(sunLonS3_deg - dipLon_deg, 360));
                
                % Fetch model coefficients
                if strcmp(MPmodel, 'Bode1994')
                    % Time-dependent coefficients
                    Gnm = coeffsBode1994G(alpha);
                else
                    % Fixed orientation coefficients
                    GnmStatic = dlmread(fullfile([coeffPath 'coeffsJupiter' MPmodel 'G.csv']), ',', nHeadLines, 0);
                    Gnm = repmat(GnmStatic, 1,1,npts);
                end
                
                % Scale normalization coefficient based from Cn = 13.5 for
                % Rss = 60RJ as assumed in Engle (1992)
                Cn = 13.5 * 60^3 ./ Rss_Rp.^3;
                
                [mpBr, mpBth, mpBphi] = deal(zeros(1,npts));
                MPcoords = PSMcoords;
                thMP = acos(xyzPSM_km(3,:) / Rp_km ./ r_Rp);
                phiMP = atan2(xyzPSM_km(2,:), xyzPSM_km(1,:));

                disp(['Evaluating ' MPmodel ' magnetopause model.'])
                for n=1:min([10, Nmax])
                    
                    % Handle m=0 separately to avoid unneccessary
                    % calculations for Bphi, which is always 0 for m=0
                    Pn0 =   LegendreS(n, 0, thMP);
                    dPn0 = dLegendreS(n, 0, thMP);
                    
                    dBr = n*(r_Rp./Rss_Rp).^(n-1) .* squeeze(Gnm(n,1,:))' .*  Pn0;
                    dBth =  (r_Rp./Rss_Rp).^(n-1) .* squeeze(Gnm(n,1,:))' .* dPn0;
                    mpBr = mpBr + dBr;
                    mpBth = mpBth + dBth;

                    for m=1:n
                        
                        % Get Schmidt semi-normalized Legendre functions
                        thMPadj = thMP;
                        thMPadj(thMPadj == 0) = eps;
                        Pnm =   LegendreS(n, m, thMPadj);
                        dPnm = dLegendreS(n, m, thMP);
                        
                        dBr =    n*(r_Rp./Rss_Rp).^(n-1) .* squeeze(Gnm(n,m+1,:))' .*  Pnm .* cos(m*phiMP);
                        dBth =     (r_Rp./Rss_Rp).^(n-1) .* squeeze(Gnm(n,m+1,:))' .* dPnm .* cos(m*phiMP);
                        dBphi = -m*(r_Rp./Rss_Rp).^(n-1) .* squeeze(Gnm(n,m+1,:))' .*  Pnm ./ sin(thMPadj) .* sin(m*phiMP);

                        mpBr = mpBr + dBr;
                        mpBth = mpBth + dBth;
                        mpBphi = mpBphi + dBphi;
                        
                    end
                end
                
                mpBr = mpBr .* Cn;
                mpBth = mpBth .* Cn;
                mpBphi = mpBphi .* Cn;
            end
            
            mpBvec = zeros(3,npts);
            [mpBvec(1,:), mpBvec(2,:), mpBvec(3,:)] = Bsph2Bxyz(mpBr, mpBth, mpBphi, thMP, phiMP);
        else
            
            error(['MPmodel "' MPmodel '" not recognized.'])
            
        end
    
        mpBvecOut = zeros(3,npts);
        disp('Rotating back to IAU coordinates for combining vectors.')
        rotMat = cspice_pxform(MPcoords, S3coords, ets);
        parfor i=1:npts
            mpBvecOut(:,i) = rotMat(:,:,i) * mpBvec(:,i);
        end
        
        if SPHOUT
            thS3 = acos(xyz_Rp(3,:) ./ r_Rp);
            phiS3 = atan2(xyz_Rp(2,:), xyz_Rp(1,:));
            [mpBvecOut(1,:), mpBvecOut(2,:), mpBvecOut(3,:)] = Bxyz2Bsph( ...
                mpBvecOut(1,:), mpBvecOut(2,:), mpBvecOut(3,:), thS3, phiS3);
        end
        
    else
        % No magnetopause field model; evaluate in/out MP anyways
        
        mpBvecOut = zeros(size(xyz_km));
        % Default subsolar point radii
        switch(parentName)
            case 'Jupiter'
                Rss_Rp = 60;
                Rp_km = 71492;
            case 'Saturn'
                Rss_Rp = 25;
                Rp_km = 60268;
            case 'Uranus'
                Rss_Rp = 60;
                Rp_km = 25559;
            case 'Neptune'
                Rss_Rp = 60;
                Rp_km = 24764;
            otherwise
                error(['parentName "' parentName '" not recognized.'])
        end
        
    end
    
    Rss_km = Rss_Rp * Rp_km;
    switch(parentName)
        case 'Jupiter'
            % Evaluate in/out of magnetopause from paraboloid of revolution
            % as described in Alexeev and Belenkaya (2005)
            OUTSIDE_MP = xyzPSM_km(3,:) > Rss_km - (xyzPSM_km(1,:).^2 + xyzPSM_km(2,:).^2)/2/Rss_km;
        case 'Saturn'
            % Use paraboloid model described in Arridge et al. (2006):
            % https://doi.org/10.1029/2005JA011574
            a = [9.1, 0.24, 0.77, -1.5];
            Rss_km = a(1) * pSW_nPa.^(-a(2)) * Rp_km;
            K = a(3) + a(4) * pSW_nPa;
            r_km = sqrt(xyz_km(1,:).^2 + xyz_km(2,:).^2 + xyz_km(3,:).^2);
            thDSZ = acos(xyzPSM_km(1,:) / r_km);
            OUTSIDE_MP = r_km > Rss_km .* (2 / (1 + cos(thDSZ))).^K;
        otherwise
            % Do the same as Jupiter
            Rss_km = Rss_Rp * Rp_km;
            OUTSIDE_MP = xyzPSM_km(3,:) > Rss_km - (xyzPSM_km(1,:).^2 + xyzPSM_km(2,:).^2)/2/Rss_km;
    end
    
end