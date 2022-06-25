moonName = 'Europa';
% Spacecraft era (sets timespan of field model)
era = 'Galileo';
CALC_NEW = 1;
ALL_MODELS = 1;
DO_FFT = 0;
outData = 'out/';

nptsApprox = 6000000;
magPhase = 0;

parentName = LoadSpice(moonName);

if CALC_NEW
    [R_P, ~, ~, ~, ~, Tparent_s, Tmoon_s] = GetBodyParams(moonName);
    Tparent_h = Tparent_s / 3600;
    Tmoon_h = Tmoon_s / 3600;
    Tsyn_h = 1/(1/Tparent_h - 1/Tmoon_h);

    if strcmp(parentName,'Saturn')
        Tinterest_h = Tmoon_h;
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
    [rM_km, latM_deg, lonM_deg] = GetPosSpice(moonName, parentName, t_h);
    altM_km = rM_km - R_P;
end

if ALL_MODELS
    switch parentName
        case 'Jupiter'; nOpts = 6;
        case 'Saturn';  nOpts = 2;
        case 'Uranus';  nOpts = 1;
        case 'Neptune'; nOpts = 1;
    end
    opts = 1:nOpts;
else
    opts = 0:0;
end
for opt=opts
    [MagModel, CsheetModel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt);
    
    if CALC_NEW
        disp(['Evaluating ' magModelDescrip ' field model.'])
        if strcmp(magModelDescrip, 'Khurana & Schwarzl 2007')
            [Bvec, ~, ~] = KSMagFldJupiter(latM_deg, lonM_deg, altM_km, t_h*3600);
        else
            [Bvec, ~, ~] = MagFldParent(parentName, latM_deg, lonM_deg, altM_km, MagModel, CsheetModel, magPhase);
        end

        disp(['Rotating field vectors into IAU frame.']);
        rotMat = cspice_pxform(['IAU_' upper(parentName)], ['IAU_' upper(moonName)], t_h*3600);
        BvecMoon = zeros(3,npts);
        parfor i=1:npts
            BvecMoon(:,i) = rotMat(:,:,i) * Bvec(:,i);
        end
        Bx = BvecMoon(1,:);
        By = BvecMoon(2,:);
        Bz = BvecMoon(3,:);
        
        save(fullfile([outData 'evalB' moonName fEnd]), 't_h', 'Bx', 'By', 'Bz');
    else
        load(fullfile([outData 'evalB' moonName fEnd]));
    end

    BD = PCA_decomposition(t_h*3600, upper(moonName), Bx, By, Bz, magModelDescrip);

    T_h = 1 ./ BD.f / 3600;
    npeaks = length(T_h);
    B0x = BD.BdxO * ones(1,npeaks);
    B0y = BD.BdyO * ones(1,npeaks);
    B0z = BD.BdzO * ones(1,npeaks);
    B1x_Re = BD.Bdxi;
    B1x_Im = BD.Bdxq;
    B1y_Re = BD.Bdyi;
    B1y_Im = BD.Bdyq;
    B1z_Re = BD.Bdzi;
    B1z_Im = BD.Bdzq;

    Tmax = 200;

    if DO_FFT
        [TfinalFFT_h, B0xFFT, B0yFFT, B0zFFT, B1xFFT, B1yFFT, B1zFFT] ...
                    = ExcitationSpectrum(moonName, nOsc, rate, Tinterest_h);
    end

    npts = length(t_h);

    omega_ph = 2*pi./T_h;
    B1x = B1x_Re + 1i*B1x_Im;
    B1y = B1y_Re + 1i*B1y_Im;
    B1z = B1z_Re + 1i*B1z_Im;
    Bmag = sqrt(abs(B1x).^2 + abs(B1y).^2 + abs(B1z).^2);

    [BxReprod, ByReprod, BzReprod] = deal(zeros(npeaks,npts));

    BeStrings = strings(npeaks,1);
    disp('Primary excitation modes:');
    BeModes = sqrt(B1x_Re.^2 + B1x_Im.^2 + B1y_Re.^2 + B1y_Im.^2 + B1z_Re.^2 + B1z_Im.^2);
    sortBeModes = sort(BeModes, 'descend');
    for i=1:npeaks
        BxReprod(i,:) = real(B1x(i) * exp(-1i * omega_ph(i) * t_h));
        ByReprod(i,:) = real(B1y(i) * exp(-1i * omega_ph(i) * t_h));
        BzReprod(i,:) = real(B1z(i) * exp(-1i * omega_ph(i) * t_h));
        indB = find(BeModes == sortBeModes(i));
        BeStrings(i) = sprintf('|B| amp: %.4e   f_Hz: %.4e   T_h: %.18f', BeModes(indB), 1/3600/T_h(indB), T_h(indB));
    end
    disp(BeStrings);

    BxTot = mean(B0x) + sum(BxReprod,1);
    ByTot = mean(B0y) + sum(ByReprod,1);
    BzTot = mean(B0z) + sum(BzReprod,1);

    figure; hold on;
    set(gcf,'Name', ['Bx full time, ' magModelDescrip]);
    plot(t_h, Bx);
    plot(t_h, BxTot);
    xlabel('Time (hr)');
    ylabel('B_x (nT)');
    legend('B_x model', 'B_x exc');
    figure; hold on;
    set(gcf,'Name', ['By full time, ' magModelDescrip]);
    plot(t_h, By);
    plot(t_h, ByTot);
    xlabel('Time (hr)');
    ylabel('B_y (nT)');
    legend('B_y model', 'B_y exc');
    figure; hold on;
    set(gcf,'Name', ['Bz full time, ' magModelDescrip]);
    plot(t_h, Bz);
    plot(t_h, BzTot);
    xlabel('Time (hr)');
    ylabel('B_z (nT)');
    legend('B_z model', 'B_z exc');

    BxD = Bx - BxTot;
    ByD = By - ByTot;
    BzD = Bz - BzTot;

    BxD1 = conj(fftshift(fft(BxD)))/npts;
    ByD1 = conj(fftshift(fft(ByD)))/npts;
    BzD1 = conj(fftshift(fft(BzD)))/npts;
    dt = t_h(2) - t_h(1);
    Fsample = 1/(dt*3600);
    dF = Fsample/npts;
    f = -Fsample/2:dF:Fsample/2-dF;
    TFT_h = 1./f ./ 3600;

    iTmax = find(TFT_h >= Tmax);
    iTmax = iTmax(end);

    Tfinal_h = TFT_h(iTmax:end);
    BxD1f = 2*BxD1(iTmax:end);
    ByD1f = 2*ByD1(iTmax:end);
    BzD1f = 2*BzD1(iTmax:end);

    figure; hold on;
    set(gcf,'Name', ['Residual FFT for ' era ' era, ' magModelDescrip]);
    plot(Tfinal_h, abs(BxD1f));
    plot(Tfinal_h, abs(ByD1f));
    plot(Tfinal_h, abs(BzD1f));
    legend('\Delta B_x', '\Delta B_y', '\Delta B_z');

    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Period (hr)');
    ylabel('FT\{B_i - B_{i,exc}\} (nT)');
    xlim([1 Tmax]);

    Bdiff = sqrt(abs(BxD1f).^2 + abs(ByD1f).^2 + abs(BzD1f).^2);
    maxBdiff = sort(Bdiff, 'descend');
    nGreatest = 10;
    disp('Highest-power remaining modes:')
    remStrings = strings(nGreatest,1);
    for i=1:nGreatest
        ind = find(Bdiff == maxBdiff(i));
        remStrings(i) = sprintf('|B| amp: %.4e   f_Hz: %.4e   T_h: %.18f', Bdiff(ind), 1/3600/Tfinal_h(ind), Tfinal_h(ind));
    end
    disp(remStrings);

    outFname = fullfile([outData 'Be1xyz_' moonName '_' era '_' fEnd '.txt']);
    header = 'period(hr),B0x(nT),B0y(nT),B0z(nT),Bex_Re(nT),Bex_Im(nT),Bey_Re(nT),Bey_Im(nT),Bez_Re(nT),Bez_Im(nT)';
    spectrumData = [T_h' B0x' B0y' B0z' real(B1x)' imag(B1x)' real(B1y)' imag(B1y)' real(B1z)' imag(B1z)'];
    dlmwrite(outFname, header, 'delimiter','');
    dlmwrite(outFname, spectrumData, 'delimiter',',', 'precision',18, '-append');

    disp(['Saved Be data to file: ' outFname])
end
