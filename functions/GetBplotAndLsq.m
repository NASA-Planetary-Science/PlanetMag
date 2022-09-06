function GetBplotAndLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
        scName, parentName, S3coords, orbStr, opt, MPopt, SEQUENTIAL, jt_h)
    npts = length(t_h);
    if length(scName) > 1
        scName = strjoin(scName, '+');
    end
    
    if ~exist('jt_h', 'var'); jt_h = []; end
    if isempty(jt_h); JUNOTOO=0; else; JUNOTOO=1; end
    
    [MagModel, CsheetModel, MPmodel, magModelDescrip, ~] = GetModelOpts(parentName, opt, MPopt);
    magPhase = 0;

    Nmax = 3;
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
    
    defName = magModelDescrip;
    scDataName = strcat(scName, ' MAG');
    commonTitle = strcat(parentName, ' model comparison vs ', scName, ' data');
    if SEQUENTIAL
        xx = 1:npts;
        xDescrip = 'Measurement index';
    elseif strcmp(parentName, 'Uranus')
        xx = t_h + 122154.0036;
        xDescrip = 'Time relative to CA (h)';
    elseif strcmp(parentName, 'Neptune')
        xx = t_h + 90752.0566;
        xDescrip = 'Time relative to CA (h)';
    else
        xx = t_h;
        xDescrip = 'Time past J2000 (h)';
    end
    figure; hold on;
    set(gcf,'Name', [char(scName) 'Br, ' orbStr ', ' magModelDescrip]);
    plot(xx, Br, 'DisplayName', defName);
    plot(xx, BrSC, 'DisplayName', scDataName);
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    legend();
    figure; hold on;
    set(gcf,'Name', [char(scName) 'Bth, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bth, 'DisplayName', defName);
    plot(xx, BthSC, 'DisplayName', scDataName);
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    legend();
    figure; hold on;
    set(gcf,'Name', [char(scName) 'Bphi, ' orbStr ', ' magModelDescrip]);
    plot(xx, Bphi, 'DisplayName', defName);
    plot(xx, BphiSC, 'DisplayName', scDataName);
    xlabel(xDescrip);
    ylabel('Vector component (nT)');
    title(commonTitle);
    legend();

    BrD = Br - BrSC;
    BthD = Bth - BthSC;
    BphiD = Bphi - BphiSC;

    figure; hold on;
    set(gcf,'Name', ['Vector comp diffs, ' orbStr ', ' magModelDescrip ' - ' char(scName) ' MAG']);
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