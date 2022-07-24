moonName = 'Triton';
% Spacecraft era (sets timespan of field model)
era = 'Voyager';
coordType = 'IAU';
CALC_NEW = 1;
ALL_MODELS = 0;
DO_FFT = 0;
DO_MPAUSE = 0;
specificModel = 0; % Set this to 0 to use default, or a number to use an opt
MPopt = 0; % As above, for magnetopause models
outData = 'out/';

nptsApprox = 600000;
magPhase = 0;

switch coordType
    case 'IAU'
        coords = ['IAU_' upper(moonName)];
    case 'SPRH'
        coords = [upper(planetName) '_SPRH'];
    otherwise
        warning(['coordType "' coordType '" not recognized Defaulting to IAU.'])
        coords = ['IAU_' upper(moonName)];
end

parentName = LoadSpice(moonName);

if CALC_NEW
    [Rp_km, ~, ~, a_AU, ~, ~, Tparent_s, Tmoon_s, nutPrecParent, nutPrecMoon] ...
        = GetBodyParams(moonName);
    Tparent_h = Tparent_s / 3600;
    Tmoon_h = Tmoon_s / 3600;
    Tsyn_h = 1/(1/Tparent_h - 1/Tmoon_h);

    if strcmp(parentName,'Saturn')
        if strcmp(moonName,'Enceladus')
            Tinterest_h = 32.927200612354675;
        else
            Tinterest_h = Tmoon_h;
        end
    else
        Tinterest_h = Tsyn_h;
    end

    switch era
        case 'Galileo'
            tStart_yr = 1995.9;
            tEnd_yr = 2003.75;
        case 'Cassini'
            tStart_yr = 2004.5;
            tEnd_yr = 2017.75;
        case 'Juno'
            tStart_yr = 2016.5;
            tEnd_yr = 2024.0;
        case 'Clipper'
            tStart_yr = 2030.0;
            tEnd_yr = 2035.0;
        case 'Voyager'
            switch parentName
                case 'Uranus'
                    tStart_yr = 1985.75;
                    tEnd_yr = 1986.25;
                case 'Neptune'
                    tStart_yr = 1989.25;
                    tEnd_yr = 1989.75;
            end
    end

    tStart_h = (2000 - tStart_yr) * 365.25*24;

    approxDur_h = (tEnd_yr - tStart_yr) * 365.25*24;
    nOsc = floor(approxDur_h / Tinterest_h);
    rate = ceil(nptsApprox/nOsc);

    dt = Tinterest_h/rate;
    duration = nOsc*Tinterest_h - dt;

    tEnd_h = tStart_h + duration;
    t_h = tStart_h:dt:tEnd_h;
    npts = length(t_h);
    disp(['Getting ' moonName ' positions for ' num2str(npts) ' pts over the ' era ' era.'])
    [rM_km, latM_deg, lonM_deg, xyz_km] = GetPosSpice(moonName, parentName, t_h);
    altM_km = rM_km - Rp_km;
end

if strcmp(coords(end-3:end), 'SPRH')
    SPHOUT = 1;
else
    SPHOUT = 0;
end

if ALL_MODELS
    switch parentName
        case 'Jupiter'; nOpts = 7;
        case 'Saturn';  nOpts = 2;
        case 'Uranus';  nOpts = 1;
        case 'Neptune'; nOpts = 1;
    end
    opts = 1:nOpts;
else
    opts = specificModel:specificModel;
end
for opt=opts
    [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt, MPopt);
    
    if CALC_NEW
        disp(['Evaluating ' magModelDescrip ' field model.'])
        if strcmp(magModelDescrip, 'Khurana & Schwarzl 2005')
            [Bvec, Mdip_nT, Odip_km] = KSMagFldJupiter(latM_deg, lonM_deg, altM_km, t_h*3600, SPHOUT);
        else
            [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, latM_deg, lonM_deg, altM_km, MagModel, CsheetModel, magPhase, SPHOUT);
        end
        if DO_MPAUSE
            
            nSW_pcc = 4 / a_AU^2 * ones(1,npts);
            vSW_kms = 400  * ones(1,npts);
            [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, t_h*3600, xyz_km, Mdip_nT, ...
                               parentName, MPmodel, SPHOUT);
            Bvec = Bvec + mpBvec;
            Bvec(:,OUTSIDE_MP) = 0;
            
        end

        disp(['Rotating field vectors into ' coords ' frame.']);
        if strcmp(coords, ['IAU_' upper(moonName)])
            BvecMoon = zeros(3,npts);
            rotMat = cspice_pxform(['IAU_' upper(parentName)], coords, t_h*3600);
            parfor i=1:npts
                BvecMoon(:,i) = rotMat(:,:,i) * Bvec(:,i);
            end
        elseif strcmp(coords, [upper(parentName) '_SPRH'])
            BvecMoon = Bvec;
        else
            if SPHOUT
                coords = [upper(parentName) '_SPRH'];
            else
                coords = ['IAU_' upper(parentName)];
            end
            warning(['Unrecognized coordinates. Defaulting to ' coords '.'])
            BvecMoon = Bvec;
        end
        
        save(fullfile([outData 'evalB' moonName fEnd]), 't_h', 'BvecMoon', 'coords');
    else
        load(fullfile([outData 'evalB' moonName fEnd]));
    end

    BD = PCA_decomposition(t_h*3600, moonName, parentName, BvecMoon, ...
        magModelDescrip, SPHOUT, 1, 1);

    T_h = 1 ./ BD.f / 3600;
    npeaks = length(T_h);
    B0vec1 = BD.Bdvec1o * ones(1,npeaks);
    B0vec2 = BD.Bdvec2o * ones(1,npeaks);
    B0vec3 = BD.Bdvec3o * ones(1,npeaks);
    B1vec1_Re = BD.Bdvec1i;
    B1vec1_Im = BD.Bdvec1q;
    B1vec2_Re = BD.Bdvec2i;
    B1vec2_Im = BD.Bdvec2q;
    B1vec3_Re = BD.Bdvec3i;
    B1vec3_Im = BD.Bdvec3q;

    Tmax = 500;

    if DO_FFT
        if CALC_NEW
            [TfinalFFT_h, B0vecFFT, B1vecFFT] ...
                        = ExcitationSpectrum(moonName, nOsc, rate, Tinterest_h, SPHOUT);
        end
        PlotSpectrum(moonName);
    end

    npts = length(t_h);

    omega_ph = 2*pi./T_h;
    B1vec1 = B1vec1_Re + 1i*B1vec1_Im;
    B1vec2 = B1vec2_Re + 1i*B1vec2_Im;
    B1vec3 = B1vec3_Re + 1i*B1vec3_Im;
    Bmag = sqrt(abs(B1vec1).^2 + abs(B1vec2).^2 + abs(B1vec3).^2);

    [Bvec1Reprod, Bvec2Reprod, Bvec3Reprod] = deal(zeros(npeaks,npts));

    BeStrings = strings(npeaks,1);
    disp('Primary excitation modes:');
    BeModes = sqrt(B1vec1_Re.^2 + B1vec1_Im.^2 + B1vec2_Re.^2 + B1vec2_Im.^2 + B1vec3_Re.^2 + B1vec3_Im.^2);
    sortBeModes = sort(BeModes, 'descend');
    for i=1:npeaks
        Bvec1Reprod(i,:) = real(B1vec1(i) * exp(-1i * omega_ph(i) * t_h));
        Bvec2Reprod(i,:) = real(B1vec2(i) * exp(-1i * omega_ph(i) * t_h));
        Bvec3Reprod(i,:) = real(B1vec3(i) * exp(-1i * omega_ph(i) * t_h));
        indB = find(BeModes == sortBeModes(i));
        BeStrings(i) = sprintf('|B| amp: %.4e   f_Hz: %.4e   T_h: %.18f', BeModes(indB), 1/3600/T_h(indB), T_h(indB));
    end
    disp(BeStrings);

    Bvec1Tot = mean(B0vec1) + sum(Bvec1Reprod,1);
    Bvec2Tot = mean(B0vec2) + sum(Bvec2Reprod,1);
    Bvec3Tot = mean(B0vec3) + sum(Bvec3Reprod,1);
    if SPHOUT
        Bv1lbl = 'B_r'; Bv2lbl = 'B_\theta'; Bv3lbl = 'B_\phi';
    else
        Bv1lbl = 'B_x'; Bv2lbl = 'B_y'; Bv3lbl = 'B_z';
    end

    figure; hold on;
    set(gcf,'Name', ['Bvec1 full time, ' magModelDescrip]);
    plot(t_h, BvecMoon(1,:));
    plot(t_h, Bvec1Tot);
    xlabel('Time (hr)');
    ylabel([Bv1lbl ' (nT)']);
    legend([Bv1lbl ' model'], [Bv1lbl ' exc']);
    figure; hold on;
    set(gcf,'Name', ['Bvec2 full time, ' magModelDescrip]);
    plot(t_h, BvecMoon(2,:));
    plot(t_h, Bvec2Tot);
    xlabel('Time (hr)');
    ylabel([Bv2lbl ' (nT)']);
    legend([Bv2lbl ' model'], [Bv2lbl ' exc']);
    figure; hold on;
    set(gcf,'Name', ['Bvec3 full time, ' magModelDescrip]);
    plot(t_h, BvecMoon(3,:));
    plot(t_h, Bvec3Tot);
    xlabel('Time (hr)');
    ylabel([Bv3lbl ' (nT)']);
    legend([Bv3lbl ' model'], [Bv3lbl ' exc']);

    Bvec1D = BvecMoon(1,:) - Bvec1Tot;
    Bvec2D = BvecMoon(2,:) - Bvec2Tot;
    Bvec3D = BvecMoon(3,:) - Bvec3Tot;

    Bvec1D1 = conj(fftshift(fft(Bvec1D)))/npts;
    Bvec2D1 = conj(fftshift(fft(Bvec2D)))/npts;
    Bvec3D1 = conj(fftshift(fft(Bvec3D)))/npts;
    dt = t_h(2) - t_h(1);
    Fsample = 1/(dt*3600);
    dF = Fsample/npts;
    f = -Fsample/2:dF:Fsample/2-dF;
    TFT_h = 1./f ./ 3600;

    iTmax = find(TFT_h >= Tmax);
    iTmax = iTmax(end);

    Tfinal_h = TFT_h(iTmax:end);
    Bvec1D1f = 2*Bvec1D1(iTmax:end);
    Bvec2D1f = 2*Bvec2D1(iTmax:end);
    Bvec3D1f = 2*Bvec3D1(iTmax:end);

    figure; hold on;
    set(gcf,'Name', ['Residual FFT for ' era ' era, ' magModelDescrip]);
    plot(Tfinal_h, abs(Bvec1D1f));
    plot(Tfinal_h, abs(Bvec2D1f));
    plot(Tfinal_h, abs(Bvec3D1f));
    legend(['\Delta ' Bv1lbl], ['\Delta ' Bv2lbl], ['\Delta ' Bv3lbl]);

    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Period (hr)');
    ylabel('FT\{B_i - B_{i,exc}\} (nT)');
    xlim([1 Tmax]);

    Bdiff = sqrt(abs(Bvec1D1f).^2 + abs(Bvec2D1f).^2 + abs(Bvec3D1f).^2);
    maxBdiff = sort(Bdiff, 'descend');
    nGreatest = 10;
    disp('Highest-power remaining modes:')
    remStrings = strings(nGreatest,1);
    for i=1:nGreatest
        ind = find(Bdiff == maxBdiff(i));
        remStrings(i) = sprintf('|B| amp: %.4e   f_Hz: %.4e   T_h: %.18f', Bdiff(ind), 1/3600/Tfinal_h(ind), Tfinal_h(ind));
    end
    disp(remStrings);

    if SPHOUT
        BeType = 'Be1sph_';
        header = 'period(hr),B0r(nT),B0th(nT),B0phi(nT),Ber_Re(nT),Ber_Im(nT),Beth_Re(nT),Beth_Im(nT),Bephi_Re(nT),Bephi_Im(nT)';
    else
        BeType = 'Be1xyz_';
        header = 'period(hr),B0x(nT),B0y(nT),B0z(nT),Bex_Re(nT),Bex_Im(nT),Bey_Re(nT),Bey_Im(nT),Bez_Re(nT),Bez_Im(nT)';
    end
    outFname = fullfile([outData BeType moonName '_' era '_' fEnd '.txt']);
    spectrumData = [T_h' B0vec1' B0vec2' B0vec3' real(B1vec1)' imag(B1vec1)' real(B1vec2)' imag(B1vec2)' real(B1vec3)' imag(B1vec3)'];
    dlmwrite(outFname, header, 'delimiter','');
    dlmwrite(outFname, spectrumData, 'delimiter',',', 'precision',18, '-append');

    disp(['Saved Be data to file: ' outFname])
end
