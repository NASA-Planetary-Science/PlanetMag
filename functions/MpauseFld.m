function mpBvecOut = MpauseFld(pSW_nPa, ets, xyz_Rp, r_Rp, Mdip_nT, parentName, MPmodel, SPHOUT)
    coeffPath = './modelCoeffs/';
    nHeadLines = 2;
    npts = length(ets);
    
    if ~strcmp(MPmodel, 'None')
        
        if strcmp(MPmodel, 'AB2005')
            % Alexeev and Belenkaya (2005) magnetopause model for Jupiter:
            % https://doi.org/10.5194/angeo-23-809-2005
            
            MPcoords = [parentName(1) 'DSZ']; % Planet-Dipole-Solar-Zenith coordinates
            % Get evaluation position in PDSZ coordinates
            xyzMP_Rp = zeros(3,npts);
            disp(['Rotating to ' MPcoords ' coordinates for ' MPmodel ' magnetopause model.'])
            rotMat = cspice_pxform(['IAU_' upper(parentName)], MPcoords, ets);
            
            % AB2005 is unclear about whether psi is relative to spin axis
            % or z axis in JDSZ coords. Test with both
            BdipDSZ = zeros(3,npts);
            
            parfor i=1:npts
                xyzMP_Rp(:,i) = rotMat(:,:,i) * xyz_Rp(:,i);
                BdipDSZ(:,i) = rotMat(:,:,i) * Mdip_nT';
            end
            B0_nT = sqrt(Mdip_nT(1)^2 + Mdip_nT(2)^2 + Mdip_nT(3)^2); % Dipole field at the magnetic equator
            Rss_Rp = 39.81 ./ pSW_nPa.^(0.23); % Distance to magnetopause at subsolar point in planetary radii
            rOutDiff_Rp = r_Rp - Rss_Rp;
            if any(rOutDiff_Rp >= 0)
                warning(['AB2005 magnetopause model is only valid for r < Rss. At least one value is' ...
                    'greater, with a maximum of ' num2str(max(rOutDiff_Rp)) ' Rp greater than Rss.'])
            end
            
            psi = acos(Mdip_nT(3) / B0_nT); % Angle of dipole tilt relative to spin axis
            psiDSZ = acos(BdipDSZ(3,:) / B0_nT); % Angle of dipole relative to Sun direction
            
            coeffs = dlmread(fullfile([coeffPath 'coeffsJupiterMpauseAB2005.csv']), ',', nHeadLines, 0);
            dPara = coeffs(:,1);
            dPerp = coeffs(:,2);
            Dtail = coeffs(:,3);
            Gtail = coeffs(:,4);
            
            [mpBr, mpBth, mpBphi] = deal(zeros(1,npts));
            thMP = acos(xyzMP_Rp(3,:) ./ r_Rp);
            phiMP = atan2(xyzMP_Rp(2,:), xyzMP_Rp(1,:));
            
            disp(['Evaluating ' MPmodel ' magnetopause model.'])
            mpB0 = B0_nT./Rss_Rp.^3;
            for n=1:6
                
                % Get unnormalized Legendre functions needed for dPara, dPerp
                Pn0 = LegendreS(n, 0, thMP, 1);
                Pn1 = LegendreS(n, 1, thMP, 1);
                dPn0 = dLegendreS(n, 0, thMP, 1);
                dPn1 = dLegendreS(n, 1, thMP, 1);
                thMPadj = thMP;
                thMPadj(thMPadj == 0) = eps;
                
                dBr = n*(r_Rp./Rss_Rp).^(n-1) .* ( dPara(n)*sin(psi)*Pn0  + dPerp(n)*cos(psi)*cos(phiMP).*Pn1 );
                dBth =  (r_Rp./Rss_Rp).^(n-1) .* ( dPara(n)*sin(psi)*dPn0 + dPerp(n)*cos(psi)*cos(phiMP).*dPn1);
                dBphi = (r_Rp./Rss_Rp).^(n-1) .* (-dPerp(n)*cos(psi)*sin(phiMP).*Pn1./sin(thMPadj));
                
                mpBr = mpBr + dBr;
                mpBth = mpBth + dBth;
                mpBphi = mpBphi + dBphi;
                
            end
            mpBr = mpBr .* mpB0;
            mpBth = mpBth .* mpB0;
            mpBphi = mpBphi .* mpB0;
            
            mpBvec = zeros(3,npts);
            [mpBvec(1,:), mpBvec(2,:), mpBvec(3,:)] = Bsph2Bxyz(mpBr, mpBth, mpBphi, thMP, phiMP);
            
        else
            
            error(['MPmodel "' MPmodel '" not recognized.'])
            
        end
    
        mpBvecOut = zeros(3,npts);
        disp('Rotating back to IAU coordinates for combining vectors.')
        rotMat = cspice_pxform(MPcoords, ['IAU_' upper(parentName)], ets);
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
        
        mpBvecOut = zeros(size(xyz_Rp));
        
    end
    
end