function GetBplotAndLsqNeptune(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
        scName, parentName, S3coords, orbStr, opt, MPopt, SEQUENTIAL, exclude)
    npts = length(t_h);
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    
    [MagModel, CsheetModel, MPmodel, magModelDescrip, ~] = GetModelOpts(parentName, opt, MPopt);
    magPhase = 0;

    Nmax = 8;
    disp(['Evaluating ' magModelDescrip ' field model with Nmax = ' num2str(Nmax) '.'])
    if strcmp(magModelDescrip, 'KS2005')
        [Bvec, Mdip_nT, Odip_km] = KSMagFldJupiter(r_km, theta, phi, ets, 1);
    else
        [Bvec, Mdip_nT, Odip_km] = MagFldParent(parentName, r_km, theta, phi, MagModel, ...
                                    CsheetModel, magPhase, 1, Nmax);
    end
    if ~strcmp(MPmodel, 'None')

        nSW_pcc = 0.14 * ones(1,npts);
        vSW_kms = 400  * ones(1,npts);
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, t_h*3600, xyz_km, ...
            Mdip_nT, Odip_km, S3coords, parentName, MPmodel, 1);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

    end

    Br = Bvec(1,:);
    Bth = Bvec(2,:);
    Bphi = Bvec(3,:);
    
    defName = 'O8, NLS relative to J2000 pole (NLS)';
    INC_PDS = 1; PDSname = 'O8 with PDS trajectory';
    INC_NM3 = 1; NM3name = 'O8, NLS with J2000 pole, Nmax=3';
    INC_PDS_NM3 = 1; PDS_NM3name = 'O8 with PDS trajectory, Nmax=3';
    INC_O8 = 1;  O8name  = 'O8, NLS with 1989 pole (NLS\_O8)';
    INC_N12_O8 = 1; N12_O8name = 'O8, IAU - 12^\circ with 1989 pole (NLS\_RADEC)';
    INC_N12 = 1; N12name = 'O8, IAU - 12^\circ (NLS\_OFFSET)';
    INC_MM3 = 1; MM3name = 'O8, MoonMag Nmax=3 with PDS trajec';
    
    if INC_PDS || INC_NM3
        BrPDS = Br;
        BthPDS = Bth;
        BphiPDS = Bphi;
        [rNLS, thetaNLS, phiNLS, ~, ~] = GetPosSpice(char(scName), parentName, t_h, 'NLS');
        [BvecNLS, ~, ~] = MagFldParent(parentName, rNLS, thetaNLS, phiNLS, MagModel, CsheetModel, magPhase, 1, Nmax);
        Br = BvecNLS(1,:);
        Bth = BvecNLS(2,:);
        Bphi = BvecNLS(3,:);
    
        if INC_NM3
            [BvecNM3, ~, ~] = MagFldParent(parentName, rNLS, thetaNLS, phiNLS, MagModel, CsheetModel, magPhase, 1, 3);
            BrNM3 = BvecNM3(1,:);
            BthNM3 = BvecNM3(2,:);
            BphiNM3 = BvecNM3(3,:);
        end
    end
    
    if INC_PDS_NM3
        [BvecPDS_NM3, ~, ~] = MagFldParent(parentName, r_km, theta, phi, MagModel, CsheetModel, magPhase, 1, 3);
        BrPDS_NM3 = BvecPDS_NM3(1,:);
        BthPDS_NM3 = BvecPDS_NM3(2,:);
        BphiPDS_NM3 = BvecPDS_NM3(3,:);
    end
    
    if INC_O8
        [rO8, thetaO8, phiO8, ~, ~] = GetPosSpice(char(scName), parentName, t_h, 'NLS_O8');
        [BvecO8, ~, ~] = MagFldParent(parentName, rO8, thetaO8, phiO8, MagModel, CsheetModel, magPhase, 1, Nmax);
        BrO8 = BvecO8(1,:);
        BthO8 = BvecO8(2,:);
        BphiO8 = BvecO8(3,:);
    end
    
    if INC_N12
        [rN12, thetaN12, phiN12, ~, ~] = GetPosSpice(char(scName), parentName, t_h, 'NLS_OFFSET');
        [BvecN12, ~, ~] = MagFldParent(parentName, rN12, thetaN12, phiN12, MagModel, CsheetModel, magPhase, 1, Nmax);
        BrN12 = BvecN12(1,:);
        BthN12 = BvecN12(2,:);
        BphiN12 = BvecN12(3,:);
    end
    
    if INC_N12_O8
        [rN12_O8, thetaN12_O8, phiN12_O8, ~, ~] = GetPosSpice(char(scName), parentName, t_h, 'NLS_RADEC');
        [BvecN12_O8, ~, ~] = MagFldParent(parentName, rN12_O8, thetaN12_O8, phiN12_O8, MagModel, CsheetModel, magPhase, 1, Nmax);
        BrN12_O8 = BvecN12_O8(1,:);
        BthN12_O8 = BvecN12_O8(2,:);
        BphiN12_O8 = BvecN12_O8(3,:);
    end
    
    if INC_MM3
        fullOrbFormatSpec = '%23s%25f%25f%25f%25f%25f%25f%25f%25f%25f%[^\n\r]';
        datFile = fullfile(['MAG/Voyager 2/VG2_comp.tab']);
        fileID = fopen(datFile,'r');
        mmData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', '', 'TextType', 'char', ...
            'EndOfLine', '\r\n', 'HeaderLines', 1);
        BrMM3 = mmData{8}';
        BthMM3 = mmData{9}';
        BphiMM3 = mmData{10}';
        BrMM3(exclude) = [];
        BthMM3(exclude) = [];
        BphiMM3(exclude) = [];
    end
    
    commonTitle = ['Neptune field model comparison, Nmax=' num2str(Nmax)];
    if SEQUENTIAL
        xx = 1:npts;
        xDescrip = 'Measurement index';
    elseif strcmp(parentName, 'Neptune')
        xx = t_h + 90752.0566;
        xDescrip = 'Time relative to CA (h)';
    else
        xx = t_h;
        xDescrip = 'Time past J2000 (h)';
    end
    
    figure; hold on;
    set(gcf,'Name', ['Br, ' orbStr ', ' magModelDescrip]);
    plot(xx, Br, 'DisplayName', defName);
    if INC_PDS; plot(xx, BrPDS, 'DisplayName', PDSname); end
    if INC_NM3; plot(xx, BrNM3, 'DisplayName', NM3name); end
    if INC_PDS_NM3; plot(xx, BrPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_O8;  plot(xx, BrO8, 'DisplayName', O8name); end
    if INC_N12_O8; plot(xx, BrN12_O8, 'DisplayName', N12_O8name); end
    if INC_N12; plot(xx, BrN12, 'DisplayName', N12name); end
    if INC_MM3; plot(xx, BrMM3, 'DisplayName', MM3name); end
    plot(xx, BrSC, 'DisplayName', 'Voyager 2 MAG');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    xlim([-1.5,1.5])
    legend();
    figure; hold on;
    set(gcf,'Name', ['Bth, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bth, 'DisplayName', defName);
    if INC_PDS; plot(xx, BthPDS, 'DisplayName', PDSname); end
    if INC_NM3; plot(xx, BthNM3, 'DisplayName', NM3name); end
    if INC_PDS_NM3; plot(xx, BthPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_O8;  plot(xx, BthO8, 'DisplayName', O8name); end
    if INC_N12_O8; plot(xx, BthN12_O8, 'DisplayName', N12_O8name); end
    if INC_N12; plot(xx, BthN12, 'DisplayName', N12name); end
    if INC_MM3; plot(xx, BthMM3, 'DisplayName', MM3name); end
    plot(xx, BthSC, 'DisplayName', 'Voyager 2 MAG');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    xlim([-1.5,1.5])
    legend();
    figure; hold on;
    set(gcf,'Name', ['Bphi, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bphi, 'DisplayName', defName);
    if INC_PDS; plot(xx, BphiPDS, 'DisplayName', PDSname); end
    if INC_PDS_NM3; plot(xx, BphiPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_NM3; plot(xx, BphiNM3, 'DisplayName', NM3name); end
    if INC_O8;  plot(xx, BphiO8, 'DisplayName', O8name); end
    if INC_N12_O8; plot(xx, BphiN12_O8, 'DisplayName', N12_O8name); end
    if INC_N12; plot(xx, BphiN12, 'DisplayName', N12name); end
    if INC_MM3; plot(xx, BphiMM3, 'DisplayName', MM3name); end
    plot(xx, BphiSC, 'DisplayName', 'Voyager 2 MAG');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    xlim([-1.5,1.5])
    legend();
    
    % Evaluate magnitudes
    Bmag = sqrt(Br.^2 + Bth.^2 + Bphi.^2);
    BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
    if INC_PDS; BmagPDS = sqrt(BrPDS.^2 + BthPDS.^2 + BphiPDS.^2); end
    if INC_PDS_NM3; BmagPDS_NM3 = sqrt(BrPDS_NM3.^2 + BthPDS_NM3.^2 + BphiPDS_NM3.^2); end
    if INC_NM3; BmagNM3 = sqrt(BrNM3.^2 + BthNM3.^2 + BphiNM3.^2); end
    if INC_O8;  BmagO8 = sqrt(BrO8.^2 + BthO8.^2 + BphiO8.^2); end
    if INC_N12_O8; BmagN12_O8 = sqrt(BrN12_O8.^2 + BthN12_O8.^2 + BphiN12_O8.^2); end
    if INC_N12; BmagN12 = sqrt(BrN12.^2 + BthN12.^2 + BphiN12.^2); end
    if INC_MM3; BmagMM3 = sqrt(BrMM3.^2 + BthMM3.^2 + BphiMM3.^2); end
    
    % Plot with same axes as past comparisons
    figure; hold on;
    set(gcf,'Name', ['Bmag, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bmag, 'DisplayName', defName);
    if INC_PDS; plot(xx, BmagPDS, 'DisplayName', PDSname); end
    if INC_PDS_NM3; plot(xx, BmagPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_NM3; plot(xx, BmagNM3, 'DisplayName', NM3name); end
    if INC_O8;  plot(xx, BmagO8, 'DisplayName', O8name); end
    if INC_N12_O8; plot(xx, BmagN12_O8, 'DisplayName', N12_O8name); end
    if INC_N12; plot(xx, BmagN12, 'DisplayName', N12name); end
    if INC_MM3; plot(xx, BmagMM3, 'DisplayName', MM3name); end
    plot(xx, BmagSC, 'DisplayName', 'Voyager 2 MAG');
    set(gca, 'YScale', 'log');
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    legend();
    grid on;
    xlim([-1.5,1.5])
    ylim([1e2,1e5])
    
    % Plot with same axes as Connerney et al. (1991)
    figure; hold on;
    set(gcf,'Name', ['Bmag (C1991), ' orbStr ', ' magModelDescrip]);
    plot(xx, Bmag, 'DisplayName', defName);
    if INC_PDS; plot(xx, BmagPDS, 'DisplayName', PDSname); end
    if INC_PDS_NM3; plot(xx, BmagPDS_NM3, 'DisplayName', PDS_NM3name); end
    if INC_NM3; plot(xx, BmagNM3, 'DisplayName', NM3name); end
    if INC_O8;  plot(xx, BmagO8, 'DisplayName', O8name); end
    if INC_N12_O8; plot(xx, BmagN12_O8, 'DisplayName', N12_O8name); end
    if INC_N12; plot(xx, BmagN12, 'DisplayName', N12name); end
    if INC_MM3; plot(xx, BmagMM3, 'DisplayName', MM3name); end
    plot(xx, BmagSC, 'DisplayName', 'Voyager 2 MAG');
    set(gca, 'YScale', 'log');
    xlabel('Time of day 1989 Aug 25');
    delt = 1 - (55*60+40.076)/3600;
    xticks([(-6:-1)/6 + delt, 0, (0:6)/6 + delt])
    xticklabels({'03:00','03:10','03:20','03:30','03:40','03:50','CA','04:00', ...
        '04:10','04:20','04:30','04:40','04:50','05:00'})
    ylabel('Vector component (nT)');
    title('Field magnitude comparisons, axes as in C1991');
    legend();
    grid on;
    xlim([-1+delt 1+delt])
    ylim([1e2,1e5])

    BrD = Br - BrSC;
    BthD = Bth - BthSC;
    BphiD = Bphi - BphiSC;
    BmagD = Bmag - BmagSC;

    figure; hold on;
    set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - MAG, ' defName]);
    plot(xx, BrD);
    plot(xx, BthD);
    plot(xx, BphiD);
    plot(xx, BmagD);
    xlabel(xDescrip);
    ylabel('Component difference (nT)');
    xlim([-1,1])
    title('Neptune field model comparison for Nmax=8, NLS relative to J2000 pole - Voyager 2 MAG')
    legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|')
    
    BrLsq = BrD.^2;
    BthLsq = BthD.^2;
    BphiLsq = BphiD.^2;
    chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3 / npts;
    disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])

    %%
    if INC_MM3
        BrD = BrNM3 - BrMM3;
        BthD = BthNM3 - BthMM3;
        BphiD = BphiNM3 - BphiMM3;
        BmagD = BmagNM3 - BmagMM3;

        figure; hold on;
        set(gcf,'Name', ['Vector comp diffs (Nmax=3), ' orbStr ', ' magModelDescrip ', PlanetMag - MoonMag']);
        plot(xx, BrD);
        plot(xx, BthD);
        plot(xx, BphiD);
        plot(xx, BmagD);
        xlabel(xDescrip);
        ylabel('Component difference (nT)');
        xlim([-1,1])
        title('Neptune field model comparison for Nmax=3, PlanetMag - MoonMag')
        legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|')
        
        BrLsq = BrD.^2;
        BthLsq = BthD.^2;
        BphiLsq = BphiD.^2;
        chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3 / npts;
        disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
    end
    %%
    if INC_N12
        BrD = BrN12 - BrSC;
        BthD = BthN12 - BthSC;
        BphiD = BphiN12 - BphiSC;
        BmagD = BmagN12 - BmagSC;

        figure; hold on;
        set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - MAG, ' N12name]);
        plot(xx, BrD);
        plot(xx, BthD);
        plot(xx, BphiD);
        plot(xx, BmagD);
        xlabel(xDescrip);
        ylabel('Component difference (nT)');
        xlim([-1,1])
        title('Neptune field model comparison, IAU with 12^\circ offset');
        legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|');
        

        BrLsq = BrD.^2;
        BthLsq = BthD.^2;
        BphiLsq = BphiD.^2;
        chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3 / npts;
        disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
    end

    %%
    if INC_O8
        BrD = BrO8 - BrSC;
        BthD = BthO8 - BthSC;
        BphiD = BphiO8 - BphiSC;
        BmagD = BmagO8 - BmagSC;

        figure; hold on;
        set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - MAG, ' N12_O8name]);
        plot(xx, BrD);
        plot(xx, BthD);
        plot(xx, BphiD);
        plot(xx, BmagD);
        xlabel(xDescrip);
        ylabel('Component difference (nT)');
        xlim([-1,1])
        title('Neptune field model comparison, NLS exactly as defined with O8 model');
        legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi', '\Delta |B|');
        
        BrLsq = BrD.^2;
        BthLsq = BthD.^2;
        BphiLsq = BphiD.^2;
        chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3 / npts;
        disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
    end
end