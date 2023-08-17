%% This doesn't work correctly yet.

function L = LshellTrace(parentName, opt, MPopt, DO_MPAUSE, sc, t_h, Nmaxin)
    if ~exist('Nmaxin', 'var'); Nmaxin = Inf; end

    [MagModel, CsheetModel, MPmodel, magModelDescrip, ~] = GetModelOpts(parentName, opt, MPopt);
    if ~DO_MPAUSE; MPmodel = 'None'; end
    [g, h, ~, ~, Rp_km, Nmax, NmaxExt] = GetGaussCoeffs(parentName, MagModel);
    Nmax = min(Nmax, Nmaxin);
    NmaxExt = min(NmaxExt, Nmaxin);
    
    ets = t_h * 3600;
    npts = length(ets);
    dipCoords = ['SM' parentName(1)];
    [~, ~, ~, xyz0_km, S3coords] = GetPosSpice(sc, parentName, t_h);
    
    res = 0.1 * Rp_km;
    L = zeros(1,npts);
    for i=1:npts
        xyzLine_km = xyz0_km(:,i);
        [~, ~, ~, xyzMag_km, ~] = GetPosSpice(sc, parentName, t_h(i), dipCoords);
        stSign = abs(xyzMag_km(3)) / xyzMag_km(3);
        forward = 1;
        while abs(xyzMag_km(3)) / xyzMag_km(3) == stSign
            r_km = sqrt(xyzLine_km(1)^2 + xyzLine_km(2)^2 + xyzLine_km(3)^2);
            theta = acos(xyzLine_km(3)/r_km);
            phi = atan2(xyzLine_km(2), xyzLine_km(1));
            if r_km < Rp_km
                break
            end
            [nextB, nextBmag] = GetB(r_km, theta, phi, xyzLine_km, ets(i), Nmax, NmaxExt, Rp_km, ...
                MagModel, CsheetModel, MPmodel, magModelDescrip, parentName, S3coords, g, h);
            if nextBmag == 0
                if xyzLine_km == xyz0_km(:,i); break; end
                forward = -1;
                xyzLine_km = xyz0_km;
                continue
            end
            xyzLine_km = xyzLine_km + forward * nextB * res / nextBmag^(1/2);
            
            rotMat = cspice_pxform(S3coords, dipCoords, ets(i));
            xyzMag_km = squeeze(pagemtimes(rotMat, xyzLine_km));
            
            
        end
        
        if nextBmag == 0
            L(i) = nan;
        elseif r_km < Rp_km
            L(i) = r_km / Rp_km / cos(asin(xyzMag_km(3) / r_km))^2;
        else
            L(i) = sqrt(xyzMag_km(1)^2 + xyzMag_km(2)^2) / Rp_km;
        end
        
    end
    

end

function [Bvec, Bmag] = GetB(r_km, theta, phi, xyz_km, ets, Nmax, NmaxExt, Rp_km, MagModel, ...
    CsheetModel, MPmodel, magModelDescrip, parentName, S3coords, g, h)
    magPhase = 0;
    
    if strcmp(magModelDescrip, 'KS2005')
        [Bvec, ~, ~] = KSMagFldJupiter(r_km, theta, phi, ets, 1);
    else
        Bvec = MagFldParentSingle(g, h, r_km, theta, phi, Rp_km, MagModel, ...
                       CsheetModel, magPhase, Nmax, NmaxExt);
    end
    if ~strcmp(MPmodel, 'None')

        nSW_pcc = 0.14 * ones(1,npts);
        vSW_kms = 400  * ones(1,npts);
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, ets, xyz_km, ...
            Mdip_nT, Odip_km, S3coords, parentName, MPmodel, 1);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

    end
    
    Bmag = sqrt(Bvec(1)^2 + Bvec(2)^2 + Bvec(3)^2);
    
end