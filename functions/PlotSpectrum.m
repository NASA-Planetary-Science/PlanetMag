function PlotSpectrum(moonName, LIVE_PLOTS, dataDir, fPattern)
% Plot a pre-saved spectrum of magnetic oscillations as derived from an FFT.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of the moon for which to plot a spectrum. An FFT data file must be present.
% LIVE_PLOTS : bool, default=0
%   Whether to load interactive figure windows for plots (true) or print them to disk (false).
% dataDir : char, 1xD, default='out'
%   Directory where data files are saved.
% fPattern : char, 1xE, default='FTdata'
%   Pattern used for saved FFT data file names.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('dataDir', 'var'); dataDir = 'out'; end
    if ~exist('fPattern', 'var'); fPattern = 'FTdata'; end
    load(fullfile(outData, [moonName 'FTdata']), 'B1vec1', 'B1vec2', 'B1vec3', 'B1mag', 'T_h', ...
         'f_Hz', 'coordType', 'SPHOUT', 'magModelDescrip', 'Tmax', 'Tinterest_h');
    LoadSpice(moonName);
    [~, ~, ~, ~, ~, ~, Tparent_s, Torb_s, ~, ~] = GetBodyParams(moonName);
    
    figure;
    hold on;
    interpreter = 'tex';
    font = 'STIX Two Math';
    set(0,'defaulttextinterpreter',interpreter);  
    set(0,'defaultAxesTickLabelInterpreter',interpreter);  
    set(0,'defaultLegendInterpreter',interpreter);
    set(0,'defaultTextFontName',font);
    set(0,'defaultAxesFontName',font);
    set(0,'defaultLegendFontName',font);
    nm   = '\rm{}';
    bnm  = '\rm\bf\fontname{STIX Two Text}';
    math  = '\it{}';
    
    plot(T_h, abs(B1vec1), 'Color', 'b');
    plot(T_h, abs(B1vec2), 'Color', 'k');
    plot(T_h, abs(B1vec3), 'Color', [0 0.8 0]);
    plot(T_h, abs(B1mag), 'Color', 'r');
    Tparent_h = Tparent_s / 3600;
    Torb_h = Torb_s / 3600;
    Tsyn_h = 1/(1/Tparent_h - 1/Torb_h);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([2 Tmax]);
    %vline(Tsyn_h,'r');
    %vline(Torb_h,'r');
    %vline(Tinterest_h,'r');
    xlabel('Period (h)');
    ylabel('Amplitude (nT)');
    set(gcf,'Name', [moonName ' FFT excitation spectrum, ' magModelDescrip]);
    title([bnm moonName ' magnetic excitation spectrum, ' magModelDescrip ', ' coordType ...
          ' coordinates'], 'fontsize', 16);
    if SPHOUT
        Bv1lbl = 'B_r'; Bv2lbl = 'B_\theta'; Bv3lbl = 'B_\phi';
    else
        Bv1lbl = 'B_x'; Bv2lbl = 'B_y'; Bv3lbl = 'B_z';
    end
    legend({[math Bv1lbl], [math Bv2lbl], [math Bv3lbl], ['|' bnm 'B' nm '|']}, 'Location', ...
           'Southeast', 'Interpreter', 'tex');
end
