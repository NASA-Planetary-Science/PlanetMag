function GetBplotAndLsq(ets, t_h, alt_km, lat_deg, lon_deg, BrSC, BthSC, BphiSC, ...
        scName, parentName, orbStr, opt, SEQUENTIAL)
    npts = length(t_h);
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    
    [MagModel, CsheetModel, magModelDescrip, ~] = GetModelOpts(parentName, opt);
    magPhase = 0;

    disp(['Evaluating ' magModelDescrip ' field model.'])
    if strcmp(magModelDescrip, 'Khurana & Schwarzl 2007')
        [Bvec, ~, ~] = KSMagFldJupiter(lat_deg, lon_deg, alt_km, ets, 1);
    else
        [Bvec, ~, ~] = MagFldParent(parentName, lat_deg, lon_deg, alt_km, MagModel, CsheetModel, magPhase, 1);
    end

    Br = Bvec(1,:);
    Bth = Bvec(2,:);
    Bphi = Bvec(3,:);
    
    if SEQUENTIAL
        xx = 1:npts;
        xDescrip = 'Measurement index';
    else
        xx = t_h;
        xDescrip = 'Time past J2000 (h)';
    end
    figure; hold on;
    set(gcf,'Name', ['Br, ' orbStr ', ' magModelDescrip]);
    plot(xx, Br);
    plot(xx, BrSC);
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    legend(magModelDescrip, strcat(scName, ' MAG'));
    figure; hold on;
    set(gcf,'Name', ['Bth, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bth);
    plot(xx, BthSC);
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    legend(magModelDescrip, strcat(scName, ' MAG'));
    figure; hold on;
    set(gcf,'Name', ['Bphi, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bphi);
    plot(xx, BphiSC);
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    legend(magModelDescrip, strcat(scName, ' MAG'));

    BrD = Br - BrSC;
    BthD = Bth - BthSC;
    BphiD = Bphi - BphiSC;

    figure; hold on;
    set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - MAG']);
    plot(xx, BrD);
    plot(xx, BthD);
    plot(xx, BphiD);
    xlabel(xDescrip);
    ylabel('Component difference (nT)');
    legend('\Delta B_r', '\Delta B_\theta', '\Delta B_\phi');

    BrLsq = BrD.^2;
    BthLsq = BthD.^2;
    BphiLsq = BphiD.^2;
    chi2 = sum([BrLsq, BthLsq, BphiLsq], 'all') / 3 / npts;
    disp(['Overall, for this range and model chi^2 = ' sprintf('%.2f', chi2) '.'])
end