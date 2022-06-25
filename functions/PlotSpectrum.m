function PlotSpectrum(moonName)
    outData = 'out/';
    load(fullfile([outData moonName 'FTdata']), 'B1x', 'B1y', 'B1z', 'T_h', 'f_Hz', ...
        'B1xS3', 'B1yS3', 'B1zS3', 'coordType', 'magModelDescrip', ...
        'Tmax', 'Tinterest_h');
    LoadSpice(moonName);
    [~, ~, ~, ~, ~, Tparent_s, Torb_s] = GetBodyParams(moonName);
    cspice_kclear;
    
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
    
    JS3 = 0;
    if JS3
        Bx = B1xS3;
        By = B1yS3;
        Bz = B1zS3;
    else
        Bx = B1x;
        By = B1y;
        Bz = B1z;
    end
    plot(T_h, abs(Bx), 'Color', 'b');
    plot(T_h, abs(By), 'Color', 'k');
    plot(T_h, abs(Bz), 'Color', [0 0.8 0]);
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
    title([bnm moonName ' magnetic excitation spectrum, ' magModelDescrip ', ' coordType ' coordinates'], 'fontsize', 16);
    legend({[math 'B_x'], [math 'B_y'], [math 'B_z']}, 'Location', 'Southeast', 'Interpreter', 'tex');
end