function L = LshellTrace(parentName, opt, MPopt, DO_MPAUSE, sc, t_h, Nmaxin, coeffPath)
% Determine the L shell of the target body at specified ephemeris times.
%
% Trace magnetic field line from the location of a target (e.g. a spacecraft) to the magnetic
% dipole equator to determine the radial distance in planetary radii to the planet's center of mass
% at that location, also known as the L shell.
%
% Note
% ----
% **This function is a work in progress and does not work correctly yet.**
%
% Parameters
% ----------
% parentName : char, 1xC
%   Parent body name recognized by SPICE. Sufficient kernels must be loaded in the pool to
%   determine the location of the target relative to this body at each ephemeris time.
% opt : int
%   Magnetospheric field model option selection. 0 selects the default model---see GetModelOpts for
%   more details.
% MPopt : int
%   Magnetopause current magnetic field model option selection. 0 selects the default model---see
%   GetModelOpts for more details.
% DO_MPAUSE : bool
%   Whether to include magnetopause currents in the magnetic field model.
% sc : char, 1xD
%   Target body for which to evaluate L shells.
% t_h : double, 1xN
%   Ephemeris times to evaluate in TDB hours relative to J2000.
% Nmaxin : int, default=Inf
%   Spherical harmonic degree to which to limit magnetic field model evaluation in tracing field
%   lines. Using the default (Inf) makes use of the highest degree both implemented in MagFldParent
%   and present in the model coefficients.
% coeffPath : char, 1xE, default='modelCoeffs'
%   Directory containing model coefficients files.
%
% Returns
% -------
% L : double, 1xN
%   L shell of each evaluation spacetime point for the target.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('Nmaxin', 'var'); Nmaxin = Inf; end
    if ~exist('coeffPath', 'var'); coeffPath = 'modelCoeffs'; end

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
            [nextB, nextBmag] = GetB(r_km, theta, phi, xyzLine_km, ets(i), Nmax, NmaxExt, ...
                Rp_km, MagModel, CsheetModel, MPmodel, magModelDescrip, parentName, S3coords, ...
                g, h);
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
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, ets, xyz_km, Mdip_nT, Odip_km, ...
            parentName, S3coords, MPmodel, coeffPath, 1);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

    end
    
    Bmag = sqrt(Bvec(1)^2 + Bvec(2)^2 + Bvec(3)^2);
    
end
