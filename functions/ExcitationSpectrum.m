function [Tpeak_h, B0vec, B1vec] = ...
    ExcitationSpectrum(moonName, nOsc, rate, Tinterest_h, DO_EPO, SPHOUT, DO_MPAUSE)
    
    parentName = LoadSpice(moonName);
    outData = 'out/';
    
    if ~exist('DO_EPO', 'var'); DO_EPO = 0; end
    if ~exist('nOsc', 'var'); nOsc = 5000; end
    if ~exist('rate', 'var'); rate = 100; end
    if ~exist('SPHOUT', 'var'); SPHOUT = 0; end
    
    magPhase = 0;
    
    [~, ~, ~, ~, ~, ~, Tparent_s, Tmoon_s, ~, ~] = GetBodyParams(moonName);
    Tparent_h = Tparent_s / 3600;
    Tmoon_h = Tmoon_s / 3600;
    Tsyn_h = 1/(1/Tparent_h - 1/Tmoon_h);
    Tmax = Tmoon_h * 1.5;
    
    if ~exist('Tinterest_h', 'var')
        EXPLORE = 1;
        Tinterest_h = Tsyn_h;
    else
        EXPLORE = 0;
    end
    
    opt = 0; % Use defaults for each planet
    if DO_MPAUSE
        MPopt = 0;
    else
        MPopt = -1;
    end
    [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt, MPopt);
    
    dt = Tinterest_h/rate;
    duration = nOsc*Tinterest_h - dt;
    tStart_h = 0;
    tEnd_h = tStart_h + duration;
    %t_h = tStart_h:dt:tEnd_h;
    t_h = 0:(1/3):(12*365.25*24);
    
    compare_V2021 = 0;
    if compare_V2021
        tstart = cspice_str2et('2018-09-01T00:00:00.000')/3600;
        tend = cspice_str2et('2028-09-01T00:00:00.000')/3600;
        t_h = linspace(tstart, tend, 2^19);
        dt = t_h(2) - t_h(1);
    end
    npts = length(t_h);
    disp(['Getting ' moonName ' positions for ' num2str(npts) ' pts.'])
    [rM_km, thetaM, phiM, xyz_km, S3coords] = GetPosSpice(moonName, parentName, t_h);
    
    disp(['Evaluating ' magModelDescrip ' field model for T = ' num2str(Tinterest_h) ' h.'])
    if contains(magModelDescrip, 'KS2005')
        [Bvec, Mdip_nT, Odip_km] = KSMagFldJupiter(rM_km, thetaM, phiM, t_h*3600, SPHOUT);
    else
        [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, rM_km, thetaM, phiM, MagModel, CsheetModel, magPhase, SPHOUT);
    end
    if DO_MPAUSE

        nSW_pcc = 0.14 * ones(1,npts);
        vSW_kms = 400  * ones(1,npts);
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, t_h*3600, xyz_km, ...
            Mdip_nT, Odip_km, parentName, MPmodel, SPHOUT);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

    end
    
    if SPHOUT
        BvecMoon = Bvec;
    else
        IAU = ['IAU_' upper(moonName)];
        if DO_EPO && strcmp(moonName,'Europa')
            EPhiO = 'EUROPAM_EUROPA_E_PHI_O';
            moonCoords = EPhiO;
            coordType = [math 'E\phi\Omega' nm];
        else
            moonCoords = IAU;
            coordType = 'IAU';
        end
        disp(['Rotating field vectors into ' moonCoords ' frame.']);
        rotMat = cspice_pxform(S3coords, moonCoords, t_h*3600);
        BvecMat(:,1,:) = Bvec;
        BvecMoon = squeeze(pagemtimes(rotMat, BvecMat));
    end
    
    disp('Taking FFTs.');
    B1vec1 = conj(fftshift(fft(BvecMoon(1,:))))/npts;
    B1vec2 = conj(fftshift(fft(BvecMoon(2,:))))/npts;
    B1vec3 = conj(fftshift(fft(BvecMoon(3,:))))/npts;
    
    Fsample = 1/(dt*3600);
    dF = Fsample/npts;
    f = -Fsample/2:dF:Fsample/2-dF;
    T_h = 1./f ./ 3600;
    
    iTmax = find(T_h >= Tmax);
    iTmax = iTmax(end);
    
    % Clip arrays to include only positive frequencies/periods
    B1vec1 = 2*B1vec1(iTmax:end);
    B1vec2 = 2*B1vec2(iTmax:end);
    B1vec3 = 2*B1vec3(iTmax:end);
    Bmag2 = abs(B1vec1).^2 + abs(B1vec2).^2 + abs(B1vec3).^2;
    B1mag = sqrt(Bmag2);
    T_h = T_h(iTmax:end);
    f_Hz = f(iTmax:end);
    
    disp('Complete! Saving data.');
    save(fullfile([outData moonName 'FTdata']), 'B1vec1', 'B1vec2', 'B1vec3', 'B1mag', 'T_h', 'f_Hz', ...
        'coordType', 'SPHOUT', 'magModelDescrip', 'Tmax', 'Tinterest_h');
    cspice_kclear;
    
    B0vec = [mean(BvecMoon(1,:)), mean(BvecMoon(2,:)), mean(BvecMoon(3,:))];
    
    checkPeak = islocalmax(Bmag2, 'MinSeparation', 5*dF, 'SamplePoints', f_Hz);
    iPeak = [];
    tol = 1e-6;
    while isempty(iPeak)
        iPeak = find(checkPeak & abs(T_h - Tinterest_h)/Tinterest_h < tol);
        tol = tol * 10;
    end
    % If multiple values are within the tolerance range, take the one
    % closest to our expected peak height of interest
    if length(iPeak) > 1
        dist = abs(T_h(iPeak) - Tinterest_h);
        iPeak = iPeak(min(dist) == dist);
    end
    
    B1vec = [B1vec1(iPeak), B1vec2(iPeak), B1vec3(iPeak)];
    
    if EXPLORE
        Tpeak_h = T_h(checkPeak & Bmag2 > 1);
        % Ensure the orbital period appears in the excitation spectrum
        if ~any(abs(Tpeak_h - Tmoon_h)/Tmoon_h < 1e-4)
            Tpeak_h = [Tmoon_h Tpeak_h];
        end
        
        save(fullfile([outData moonName 'TseriesData' fEnd]), 'BvecMoon', 't_h', 'Tmax');
    else
        Tpeak_h = T_h(iPeak);
    end
end