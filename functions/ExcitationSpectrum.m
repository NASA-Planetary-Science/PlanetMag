function [Tpeak_h, B0x, B0y, B0z, B1x_peak, B1y_peak, B1z_peak] = ...
    ExcitationSpectrum(moonName, nOsc, rate, Tinterest_h, DO_EPO)
    
    parentName = LoadSpice(moonName);
    outData = 'out/';
    
    if ~exist('DO_EPO', 'var'); DO_EPO = 0; end
    if ~exist('nOsc', 'var'); nOsc = 5000; end
    if ~exist('rate', 'var'); rate = 100; end
    
    magPhase = 0;
    
    [R_P, ~, ~, ~, ~, Tparent_s, Tmoon_s] = GetBodyParams(moonName);
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
    
    switch parentName
        case 'Jupiter'
            opt = 5;
        case 'Saturn'
            opt = 2;
        case 'Uranus'
            opt = 1;
        case 'Neptune'
            opt = 1;
    end
    [MagModel, CsheetModel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt);
    
    dt = Tinterest_h/rate;
    duration = nOsc*Tinterest_h - dt;
    tStart_h = 0;
    tEnd_h = tStart_h + duration;
    t_h = tStart_h:dt:tEnd_h;
    
    compare_V2021 = 0;
    if compare_V2021
        tstart = cspice_str2et('2018-09-01T00:00:00.000')/3600;
        tend = cspice_str2et('2028-09-01T00:00:00.000')/3600;
        t_h = linspace(tstart, tend, 2^19);
        dt = t_h(2) - t_h(1);
    end
    npts = length(t_h);
    disp(['Getting ' moonName ' positions for ' num2str(npts) ' pts.'])
    [rM_km, latM_deg, lonM_deg] = GetPosSpice(moonName, parentName, t_h);
    altM_km = rM_km - R_P;
    
    disp(['Evaluating ' magModelDescrip ' field model for T = ' num2str(Tinterest_h) ' h.'])
    if strcmp(magModelDescrip, 'Khurana & Schwarzl 2007')
        [Bvec, ~, ~] = kkMagFldJupiter(latM_deg, lonM_deg, altM_km, t_h*3600);
    else
        [Bvec, ~, ~] = MagFldParent(parentName, latM_deg, lonM_deg, altM_km, MagModel, CsheetModel, magPhase);
    end
    
    % Separate out components
    BxS3 = Bvec(1,:);
    ByS3 = Bvec(2,:);
    BzS3 = Bvec(3,:);
        
    spkS3 = ['IAU_' upper(parentName)];
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
    rotMat = cspice_pxform(spkS3, moonCoords, t_h*3600);
    BvecMoon = zeros(3,npts);
    parfor i=1:npts
        BvecMoon(:,i) = rotMat(:,:,i) * Bvec(:,i);
    end
    Bx = BvecMoon(1,:);
    By = BvecMoon(2,:);
    Bz = BvecMoon(3,:);
    
    disp('Taking FFTs.');
    B1x = conj(fftshift(fft(Bx)))/npts;
    B1y = conj(fftshift(fft(By)))/npts;
    B1z = conj(fftshift(fft(Bz)))/npts;
    
    B1xS3 = conj(fftshift(fft(BxS3)))/npts;
    B1yS3 = conj(fftshift(fft(ByS3)))/npts;
    B1zS3 = conj(fftshift(fft(BzS3)))/npts;
    
    Fsample = 1/(dt*3600);
    dF = Fsample/npts;
    f = -Fsample/2:dF:Fsample/2-dF;
    T_h = 1./f ./ 3600;
    
    iTmax = find(T_h >= Tmax);
    iTmax = iTmax(end);
    
    % Clip arrays to include only positive frequencies/periods
    B1x = 2*B1x(iTmax:end);
    B1y = 2*B1y(iTmax:end);
    B1z = 2*B1z(iTmax:end);
    B1xS3 = 2*B1xS3(iTmax:end);
    B1yS3 = 2*B1yS3(iTmax:end);
    B1zS3 = 2*B1zS3(iTmax:end);
    T_h = T_h(iTmax:end);
    f_Hz = f(iTmax:end);
    
    disp('Complete! Saving data.');
    save(fullfile([outData moonName 'FTdata']), 'B1x', 'B1y', 'B1z', 'T_h', 'f_Hz', ...
        'B1xS3', 'B1yS3', 'B1zS3', 'coordType', 'magModelDescrip', ...
        'Tmax', 'Tinterest_h');
    cspice_kclear;
    
    B0x = mean(Bx);
    B0y = mean(By);
    B0z = mean(Bz);
    
    Bmag2 = abs(B1x).^2 + abs(B1y).^2 + abs(B1z).^2;
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
    B1x_peak = B1x(iPeak);
    B1y_peak = B1y(iPeak);
    B1z_peak = B1z(iPeak);
    
    if EXPLORE
        Tpeak_h = T_h(checkPeak & Bmag2 > 1);
        % Ensure the orbital period appears in the excitation spectrum
        if ~any(abs(Tpeak_h - Tmoon_h)/Tmoon_h < 1e-4)
            Tpeak_h = [Tmoon_h Tpeak_h];
        end
        
        save(fullfile([outData moonName 'TseriesData' fEnd]), 'Bx', 'By', 'Bz', 'BxS3', 'ByS3', ...
            'BzS3', 't_h', 'Tmax');
    else
        Tpeak_h = T_h(iPeak);
    end
end