function [T_h, B0vec, B1vec1, B1vec2, B1vec3, outFname, header] = PlanetMag(moonName, era, ...
    coordType, CALC_NEW, ALL_MODELS, DO_FFT, DO_MPAUSE, LIVE_PLOTS, nptsApprox, specificModel, ...
    specificMPmodel, outData, figDir, coeffPath, fPatternFT, fPatternTseries, figXtn, magPhase_deg)
% Evaluates planetary magnetic field for a time series at the specified moon location and inverts
% for the complex amplitudes of oscillation in that moon's frame.
%
% This function is intended to be run as a script, so all parameters are optional. 
%
% Parameters
% ----------
% moonName : char, 1xC, default='Europa'
%   Name of target moon for which to generate magnetic spectrum amplitudes.
% era : char, 1xD, default='Galileo'
%   Time period over which measurements will be evaluated, corresponding to the prime mission
%   lifetime. Options:
%
%     - ``'Swarm'``
%     - ``'Galileo'``
%     - ``'Cassini'``
%     - ``'Juno'``
%     - ``'Clipper'``
%     - ``'Voyager'``
%
% coordType : char, 1xE, default='IAU'
%   Desired standard coordinates for magnetic spectrum inversion. Options:
% 
%     - ``'IAU'``
%     - ``'SPRH'``
%
% CALC_NEW : bool, default=1
%   Whether to perform calculations or attempt to reload saved data for plotting purposes.
% ALL_MODELS : bool, default=0
%   Whether to run all implemented model options for the desired body.
% DO_FFT : bool, default=0
%   Whether to calculate and print an FFT from the time series after performing the inversion.
% DO_MPAUSE : bool, default=0
%   Whether to include a magnetopause screening current model.
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% nptsApprox : int, default=12*365.25*3*24
%   Desired number of points to use in time series for inversion. A whole number of the period of 
%   interest (typically synodic period, as it is the strongest oscillation) will ultimately be 
%   selected, which is why this number is approximate.
% specificModel : int, default=0
%   Index number for the magnetospheric model to run. Options depend on the body, and setting to
%   0 selects the default model. See GetModelOpts for a description of the options.
% specificMPmodel : int, default=0
%   Index number for the magnetopause model to run if DO_MPAUSE is true. 0 selects the default
%   model See MpauseFld for a description of the options.
% outData : char, 1xF, default='out'
%   Directory to use for output of complex spectrum amplitudes.
% figDir : char, 1xG, default='figures'
%   Directory to use for output figures.
% coeffPath : char, 1xH, default='modelCoeffs'
%   Directory containing model coefficients files.
% fPatternFT : char, 1xI, default='FTdata'
%   Pattern for file names of FFT spectrum data saved to disk.
% fPatternTseries : char, 1xJ, default='TseriesData'
%   Pattern for file names of time series data saved to disk.
% figXtn : char, 1xK, default='pdf'
%   Extension to use for figures, which determines the file type.
% magPhase_deg : double, default=0
%   Arbitrary offset in degrees by which to rotate the magnetospheric field evaluation.
%
% Returns
% -------
% T_h : double, 1xP
%   Oscillation periods of significant peaks in the excitation spectrum in hours.
% B0vec : double, 1x3
%   Time-independent background magnetic field vector from evaluated time series in nT in the
%   selected coordinates.
% B1vec1, B1vec2, B1vec3 : complex double, 1xP
%   Degree-1 complex excitation moments inverted from the evaluated time series in nT at the body
%   center. Components 1, 2, and 3 are aligned with, respectively, x, y, z or r, theta, phi
%   depending on the selected coordinates.
% outFname : char, 1xH
%   Output file name to which the excitations moments were written.
% header : char, 1xI
%   Column header printed to output file, which contains information about what data were saved.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('moonName', 'var'); moonName = 'Europa'; end
    if ~exist('era', 'var'); era = 'Galileo'; end
    if ~exist('coordType', 'var'); coordType = 'IAU'; end
    if ~exist('CALC_NEW', 'var'); CALC_NEW = 1; end
    if ~exist('ALL_MODELS', 'var'); ALL_MODELS = 0; end
    if ~exist('DO_FFT', 'var'); DO_FFT = 0; end
    if ~exist('DO_MPAUSE', 'var'); DO_MPAUSE = 0; end
    if ~exist('LIVE_PLOTS', 'var'); LIVE_PLOTS = 0; end
    if ~exist('nptsApprox', 'var'); nptsApprox = 12*365.25*3*24; end
    if ~exist('specificModel', 'var'); specificModel = 0; end
    if ~exist('specificMPmodel', 'var'); specificMPmodel = 0; end
    if ~exist('outData', 'var'); outData = 'out'; end
    if ~exist('figDir', 'var'); figDir = 'figures'; end
    if ~exist('coeffpath', 'var'); coeffPath = 'modelCoeffs'; end
    if ~exist('fPatternFT', 'var'); fPatternFT = 'FTdata'; end
    if ~exist('fPatternTseries', 'var'); fPatternTseries = 'TseriesData'; end
    if ~exist('figXtn', 'var'); figXtn = 'pdf'; end
    if ~exist('magPhase_deg', 'var'); magPhase_deg = 0; end
    
    parentName = LoadSpice(moonName);
    SetPlotDefaults();
    % The following are defined in SetPlotDefaults. Do NOT reset them anywhere else.
    global nmTxt
    global bnmTxt
    global mathTxt
    global bmathTxt
    
    switch coordType
        case 'IAU'
            outCoords = ['IAU_' upper(moonName)];
        case 'SPRH'
            outCoords = [upper(parentName) '_SPRH'];
        otherwise
            warning(['coordType "' coordType '" not recognized. Defaulting to IAU.'])
            outCoords = ['IAU_' upper(moonName)];
    end
    
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
            case 'Swarm'
                tStart_yr = 2020.0;
                tEnd_yr = 2025.0;
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
                    case 'Jupiter'
                        tStart_yr = 1979.0;
                        tEnd_yr = 1980.0;
                    case 'Saturn'
                        tStart_yr = 1980.5;
                        tEnd_yr = 1982.0;
                    case 'Uranus'
                        tStart_yr = 1985.75;
                        tEnd_yr = 1986.25;
                    case 'Neptune'
                        tStart_yr = 1989.25;
                        tEnd_yr = 1989.75;
                end
            otherwise
                if ~strcmp(era, 'None')
                    warning(['era "' era '" not recognized. Starting at J2000 with 20 min steps.'])
                end
                tStart_yr = 2000;
                tEnd_yr = tStart_yr + nptsApprox / (365.25*24*3);
        end
    
        tStart_h = (tStart_yr - 2000) * 365.25*24;
    
        approxDur_h = (tEnd_yr - tStart_yr) * 365.25*24;
        nOsc = floor(approxDur_h / Tinterest_h);
        rate = ceil(nptsApprox/nOsc);
    
        dt = Tinterest_h/rate;
        duration = nOsc*Tinterest_h - dt;
    
        tEnd_h = tStart_h + duration;
        t_h = tStart_h:dt:tEnd_h;
        npts = length(t_h);
        disp(['Getting ' moonName ' positions for ' num2str(npts) ' pts over the ' era ' era.'])
        [rM_km, thetaM, phiM, xyz_km, S3coords] = GetPosSpice(moonName, parentName, t_h);
    end
    
    if strcmp(outCoords(end-3:end), 'SPRH')
        SPHOUT = 1;
    else
        SPHOUT = 0;
    end
    
    if ALL_MODELS
        switch parentName
            case 'Earth'; nOpts = 1; nMPopts = 0;
            case 'Jupiter'; nOpts = 7; nMPopts = 5;
            case 'Saturn';  nOpts = 2; nMPopts = 0;
            case 'Uranus';  nOpts = 2; nMPopts = 2;
            case 'Neptune'; nOpts = 1; nMPopts = 1;
        end
        opts = 1:nOpts;
        MPopts = 1:(nMPopts + 1); % Add 1 to force noMP model in addition
    else
        opts = specificModel:specificModel;
        if DO_MPAUSE
            MPopts = specificMPmodel:specificMPmodel;
        else
            MPopts = -1:-1;
        end
    end
    for opt=opts
        for MPopt=MPopts
            [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, ...
                opt, MPopt);
    
            if CALC_NEW
                disp(['Evaluating ' magModelDescrip ' field model.'])
                if contains(magModelDescrip, 'KS2005')
                    if MPopt == MPopts(1)
                        [Bvec, Mdip_nT, Odip_km] = MagFldJupiterKS2005(rM_km, thetaM, phiM, ...
                            t_h*3600, SPHOUT);
                    else
                        continue
                    end
                else
                    [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, rM_km, thetaM, phiM, ...
                        MagModel, CsheetModel, magPhase_deg, SPHOUT);
                end
                if DO_MPAUSE && ~(strcmp(MPmodel, 'None') || contains(magModelDescrip, 'KS2005'))
    
                    nSW_pcc = 4 / a_AU^2 * ones(1,npts);
                    vSW_kms = 400  * ones(1,npts);
                    [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, t_h*3600, xyz_km, ...
                        Mdip_nT, Odip_km, parentName, S3coords, MPmodel, coeffPath, SPHOUT);
                    Bvec = Bvec + mpBvec;
                    Bvec(:,OUTSIDE_MP) = 0;
    
                end
    
                disp(['Rotating field vectors into ' outCoords ' frame.']);
                if strcmp(outCoords, ['IAU_' upper(moonName)])
                    rotMat = cspice_pxform(S3coords, outCoords, t_h*3600);
                    BvecMat = zeros(3,1,npts);
                    BvecMat(:,1,:) = Bvec;
                    BvecMoon = squeeze(pagemtimes(rotMat, BvecMat));
    
                elseif strcmp(outCoords, [upper(parentName) '_SPRH'])
                    BvecMoon = Bvec;
    
                else
                    if SPHOUT
                        outCoords = [upper(parentName) '_SPRH'];
                    else
                        outCoords = S3coords;
                    end
                    warning(['Unrecognized coordinates. Defaulting to ' outCoords '.'])
                    BvecMoon = Bvec;
                end
    
                save(fullfile(outData, ['evalB' moonName fEnd]), 't_h', 'BvecMoon', 'outCoords');
            else
                load(fullfile(outData, ['evalB' moonName fEnd]));
            end
    
            BD = ICAdecomposition(moonName, parentName, t_h*3600, BvecMoon, magModelDescrip, ...
                SPHOUT, 1, 1, 1, LIVE_PLOTS, figDir, figXtn);
    
            T_h = 1 ./ BD.f / 3600;
            npeaks = length(T_h);
            B0vec1 = BD.BexcVec1o * ones(1,npeaks);
            B0vec2 = BD.BexcVec2o * ones(1,npeaks);
            B0vec3 = BD.BexcVec3o * ones(1,npeaks);
            B0vec = [BD.BexcVec1o, BD.BexcVec2o, BD.BexcVec3o];
            B1vec1_Re = BD.BexcVec1i;
            B1vec1_Im = BD.BexcVec1q;
            B1vec2_Re = BD.BexcVec2i;
            B1vec2_Im = BD.BexcVec2q;
            B1vec3_Re = BD.BexcVec3i;
            B1vec3_Im = BD.BexcVec3q;
    
            Tmax = 500;
    
            if DO_FFT
                if ~CALC_NEW
                    [Rp_km, ~, ~, a_AU, ~, ~, Tparent_s, Tmoon_s, nutPrecParent, nutPrecMoon] ...
                        = GetBodyParams(moonName);
                    Tparent_h = Tparent_s / 3600;
                    Tmoon_h = Tmoon_s / 3600;
                    Tsyn_h = 1/(1/Tparent_h - 1/Tmoon_h);
    
                    if strcmp(parentName, 'Saturn')
                        if strcmp(moonName, 'Enceladus')
                            Tinterest_h = 32.927200612354675;
                        else
                            Tinterest_h = Tmoon_h;
                        end
                    else
                        Tinterest_h = Tsyn_h;
                    end
                    tStart_h = (tStart_yr - 2000) * 365.25*24;
                    approxDur_h = (tEnd_yr - tStart_yr) * 365.25*24;
                    nOsc = floor(approxDur_h / Tinterest_h);
                end
                [TfinalFFT_h, B0vecFFT, B1vecFFT] = ExcitationSpectrum(moonName, nOsc, rate, ...
                    Tinterest_h, coordType, outData, coeffPath, fPatternFT, fPatternTseries, ...
                    magPhase_deg, 0, SPHOUT, DO_MPAUSE);
                figNum = 200 + 10*specificModel + 1*specificMPmodel;
                PlotSpectrum(moonName, LIVE_PLOTS, figNum, outData, figDir, fPatternFT, figXtn);
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
            BeModes = sqrt(B1vec1_Re.^2 + B1vec1_Im.^2 + B1vec2_Re.^2 + B1vec2_Im.^2 + ...
                B1vec3_Re.^2 + B1vec3_Im.^2);
            sortBeModes = sort(BeModes, 'descend');
            for i=1:npeaks
                Bvec1Reprod(i,:) = real(B1vec1(i) * exp(-1i * omega_ph(i) * t_h));
                Bvec2Reprod(i,:) = real(B1vec2(i) * exp(-1i * omega_ph(i) * t_h));
                Bvec3Reprod(i,:) = real(B1vec3(i) * exp(-1i * omega_ph(i) * t_h));
                indB = find(BeModes == sortBeModes(i));
                BeStrings(i) = sprintf('|B| amp: %.4e   f_Hz: %.4e   T_h: %.18f', ...
                    BeModes(indB), 1/3600/T_h(indB), T_h(indB));
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
            figNumBase = 1000 + 100*specificModel + 10*specificMPmodel;
    
            % Plot vector component 1 reproduction vs input time series
            xx = t_h;
            yy = [BvecMoon(1,:); Bvec1Tot];
            windowName = ['Bvec1 full time, ' magModelDescrip];
            legendStrings = [string([mathTxt Bv1lbl nmTxt ' model']), string([mathTxt Bv1lbl ...
                nmTxt ' exc'])];
            titleInfo = [mathTxt Bv1lbl bnmTxt ' reproduced vs. ' magModelDescrip];
            xInfo = [nmTxt 'Time (hr)'];
            yInfo = [mathTxt Bv1lbl nmTxt ' (nT)'];
            fName = [moonName 'Bv1_fullRepro' magModelDescrip era coordType];
            PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, ...
                fName, figDir, figXtn, LIVE_PLOTS, figNumBase + 1);

            % Plot vector component 2 reproduction vs input time series
            yy = [BvecMoon(2,:); Bvec2Tot];
            windowName = ['Bvec2 full time, ' magModelDescrip];
            legendStrings = [string([mathTxt Bv2lbl nmTxt ' model']), string([mathTxt Bv2lbl ...
                nmTxt ' exc'])];
            titleInfo = [mathTxt Bv2lbl bnmTxt ' reproduced vs. ' magModelDescrip];
            yInfo = [mathTxt Bv2lbl nmTxt ' (nT)'];
            fName = [moonName 'Bv2_fullRepro' magModelDescrip era coordType];
            PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, ...
                fName, figDir, figXtn, LIVE_PLOTS, figNumBase + 2);

            % Plot vector component 3 reproduction vs input time series
            yy = [BvecMoon(3,:); Bvec3Tot];
            windowName = ['Bvec3 full time, ' magModelDescrip];
            legendStrings = [string([mathTxt Bv3lbl nmTxt ' model']), string([mathTxt Bv3lbl ...
                nmTxt ' exc'])];
            titleInfo = [mathTxt Bv3lbl bnmTxt ' reproduced vs. ' magModelDescrip];
            yInfo = [mathTxt Bv3lbl nmTxt ' (nT)'];
            fName = [moonName 'Bv3_fullRepro' magModelDescrip era coordType];
            PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, ...
                fName, figDir, figXtn, LIVE_PLOTS, figNumBase + 3);
    
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
    
            % Plot residual FFT spectrum
            xx = Tfinal_h;
            yy = [abs(Bvec1D1f); abs(Bvec2D1f); abs(Bvec3D1f)];
            windowName = ['Residual FFT for ' era ' era, ' magModelDescrip];
            legendStrings = [string([mathTxt '\Delta ' Bv1lbl]), string([mathTxt '\Delta ' ...
                Bv2lbl]), string([mathTxt '\Delta ' Bv3lbl])];
            xlims = [1 Tmax];
            xInfo = [nmTxt 'Period (hr)'];
            yInfo = [nmTxt 'FT' mathTxt '\{B_i - B_{i,\rm{exc}}\}' nmTxt ' (nT)'];
            titleInfo = [bnmTxt 'Residual FFT for ' magModelDescrip ' reproduction'];
            fName = [moonName 'ResidualFFT' magModelDescrip era coordType];
            figNum = figNumBase / 10;
            PlotGeneric(xx, yy, legendStrings, windowName, titleInfo, xInfo, yInfo, ...
                fName, figDir, figXtn, LIVE_PLOTS, figNum, 'log', 'log', xlims);
    
            % Locate main frequencies of differences in the reproduced time series
            Bdiff = sqrt(abs(Bvec1D1f).^2 + abs(Bvec2D1f).^2 + abs(Bvec3D1f).^2);
            maxBdiff = sort(Bdiff, 'descend');
            nGreatest = 10;
            disp('Highest-power remaining modes:')
            remStrings = strings(nGreatest,1);
            for i=1:nGreatest
                ind = find(Bdiff == maxBdiff(i));
                remStrings(i) = sprintf('|B| amp: %.4e   f_Hz: %.4e   T_h: %.18f', Bdiff(ind), ...
                    1/3600/Tfinal_h(ind), Tfinal_h(ind));
            end
            disp(remStrings);
    
            if SPHOUT
                BeType = 'Be1sph_';
                header = ['period(hr),B0r(nT),B0th(nT),B0phi(nT),Ber_Re(nT),Ber_Im(nT),' ...
                    'Beth_Re(nT),Beth_Im(nT),Bephi_Re(nT),Bephi_Im(nT)'];
            else
                BeType = 'Be1xyz_';
                header = ['period(hr),B0x(nT),B0y(nT),B0z(nT),Bex_Re(nT),Bex_Im(nT),' ...
                    'Bey_Re(nT),Bey_Im(nT),Bez_Re(nT),Bez_Im(nT)'];
            end
            outFname = fullfile(outData, [BeType moonName '_' era '_' fEnd '.txt']);
            spectrumData = [T_h' B0vec1' B0vec2' B0vec3' real(B1vec1)' imag(B1vec1)' ...
                real(B1vec2)' imag(B1vec2)' real(B1vec3)' imag(B1vec3)'];
            dlmwrite(outFname, header, 'delimiter','');
            dlmwrite(outFname, spectrumData, 'delimiter',',', 'precision',18, '-append');
    
            disp(['Saved Be data to file: ' outFname])
        end
    end

end
