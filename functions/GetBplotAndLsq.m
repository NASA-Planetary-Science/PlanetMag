function GetBplotAndLsq(ets, t_h, alt_km, lat_deg, lon_deg, BrSC, BthSC, BphiSC, ...
        scName, parentName, orbStr, opt, MPopt, SEQUENTIAL, jt_h)
    npts = length(t_h);
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    
    if ~exist('jt_h', 'var'); jt_h = []; end
    if isempty(jt_h); JUNOTOO=0; else; JUNOTOO=1; end
    
    [MagModel, CsheetModel, MPmodel, magModelDescrip, ~] = GetModelOpts(parentName, opt, MPopt);
    magPhase = 0;

    disp(['Evaluating ' magModelDescrip ' field model.'])
    if strcmp(magModelDescrip, 'KS2005')
        [Bvec, Mdip_nT, ~] = KSMagFldJupiter(lat_deg, lon_deg, alt_km, ets, 1);
    else
        [Bvec, Mdip_nT, ~] = MagFldParent(parentName, lat_deg, lon_deg, alt_km, MagModel, ...
                                    CsheetModel, magPhase, 1);
    end
    if DO_MPAUSE

        nSW_pcc = 0.14 * ones(1,npts);
        vSW_kms = 400  * ones(1,npts);
        [mpBvec, OUTSIDE_MP] = MpauseFld(nSW_pcc, vSW_kms, t_h*3600, xyz_km, Mdip_nT, ...
                           parentName, MPmodel, SPHOUT);
        Bvec = Bvec + mpBvec;
        Bvec(:,OUTSIDE_MP) = 0;

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