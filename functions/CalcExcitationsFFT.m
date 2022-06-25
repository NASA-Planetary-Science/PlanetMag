%% NOTE: This script evaluates excitation moments with an FFT method.
%        This method is not recommended, as fine-tuning the sampling of the
%        time series data to well-resolve one peak in the Fourier spectrum
%        will necessarily contribute to spectral leakage for all other
%        non-harmonic peaks. This forces the phases of each excitation
%        moment relative to the reference time (J2000) to be incorrect.
%        Instead, the method applied using principle component
%        decomposition in PlanetMag is recommended.
%%
% Procedure:
% 0. Download SPICE kernels listed in loadSpice.
% 1. Run a spectrum for a lot of eval points, ~1M+. Since Tinterest_h is
% not set, the synodic period is used and the periods with peaks above 1 nT
% are all returned.
% 2. Run one spectrum per period of interest and save the resulting peak
% heights in x, y, and z to a text file.

moonName = 'Ganymede';
outData = 'out/';
%[Tpeaks_h, ~, ~, ~, ~, ~, ~] = ExcitationSpectrum(moonName, 5000, 200);
if strcmp(moonName, 'Europa')
    fOrb_ph = 1/(3.551181041*24);
    fSyn_ph = 1/9.92496 - fOrb_ph;
    fAnom_ph = 3.281745587e-6*3600;
    f_ph = [fSyn_ph*3, fSyn_ph*2, fSyn_ph - fOrb_ph, fSyn_ph, fSyn_ph + fOrb_ph, fAnom_ph, fOrb_ph];
elseif strcmp(moonName, 'Ganymede')
    fOrb_ph = 1/(7.154553119547032*24);
    fSyn_ph = 870.536 / 360 / 24 - fOrb_ph;
    f_ph = [fSyn_ph*3, fSyn_ph*2, fSyn_ph - fOrb_ph, fSyn_ph, fSyn_ph + fOrb_ph, fOrb_ph];
end
Tpeaks_h = 1./f_ph;
Tpeaks_h = sort(Tpeaks_h);
%plotSpectrum(moonName);
nPeaks = length(Tpeaks_h);
disp(['Found ' num2str(nPeaks) ' peaks larger than 1 nT to evaluate:' newline num2str(Tpeaks_h)])
[Tfinal_h, B0x, B0y, B0z, B1x, B1y, B1z] = deal(zeros(1,nPeaks));

% For long periods of oscillation, SPICE kernels may not contain a long 
% enough time series. For the cspice_spkpos error "Insufficient ephemeris
% data has been loaded to compute the position", tStart_h may need to be 
% reduced in ExcitationSpectrum, or nOsc can be reduced. The sampling rate
% should be increased if nOsc is decreased.
EXPLORE = 0;
if EXPLORE
    nOsc = 2000;
    rate = 3000;
    [Tfinal_h, B0x, B0y, B0z, B1x, B1y, B1z] ...
            = ExcitationSpectrum(moonName, nOsc, rate);
else
    nOsc = 4000;
    rate = 150;
    parfor i=1:nPeaks
        [Tfinal_h(i), B0x(i), B0y(i), B0z(i), B1x(i), B1y(i), B1z(i)] ...
            = ExcitationSpectrum(moonName, nOsc, rate, Tpeaks_h(i));
    end
end

if ~EXPLORE
    header = 'period(hr),B0x(nT),B0y(nT),B0z(nT),Bex_Re(nT),Bex_Im(nT),Bey_Re(nT),Bey_Im(nT),Bez_Re(nT),Bez_Im(nT)';
    spectrumData = [Tfinal_h' B0x' B0y' B0z' real(B1x)' imag(B1x)' real(B1y)' imag(B1y)' real(B1z)' imag(B1z)'];
    outFname = fullfile([outData 'Be1xyz_' moonName '.txt']);
    dlmwrite(outFname, header, 'delimiter','');
    dlmwrite(outFname, spectrumData, 'delimiter',',', 'precision',18, '-append');

    disp(['Saved Be data to file: ' outFname])
end
