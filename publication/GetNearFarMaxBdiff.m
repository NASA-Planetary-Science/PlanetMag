function GetNearFarMaxBdiff
    
    disp('Maximum near/far vector difference magnitudes are:')
    moons = [
        "Io"
        "Europa"
        "Ganymede"
        "Enceladus"
        "Mimas"
        "Miranda"
        "Ariel"
        "Triton"
        ];
    
    % Iterate over all moons
    for i=1:length(moons)
        moon = char(moons(i));
        spkMoon = upper(['IAU_' moon]);

        % Evaluate once per minute over 3 months past J2000
        ets = 0:60:86400*90;
        planet = LoadSpice(moon,'None');
        [~, ~, ~, xyz_km, spkParent] = GetPosSpice(moon, planet, ets/3600);
        [~, ~, RmoonEq_km, ~, ~, ~, ~, ~, ~, ~] = GetBodyParams(moon);

        % Transform from IAU frame (0 N, 0 E) and (0 N, 180 E) points to System III frame
        [DxNear_km, DyNear_km, DzNear_km] = RotateVecSpice(RmoonEq_km, 0, 0, ets, spkMoon, ...
            spkParent);
        [DxFar_km, DyFar_km, DzFar_km] = RotateVecSpice(-RmoonEq_km, 0, 0, ets, spkMoon, ...
            spkParent);
        DxyzNear_km = [DxNear_km; DyNear_km; DzNear_km];
        DxyzFar_km = [DxFar_km; DyFar_km; DzFar_km];
        xyzNear_km = xyz_km + DxyzNear_km;
        xyzFar_km = xyz_km + DxyzFar_km;
        [rNear_km, thetaNear, phiNear] = xyz2sph(xyzNear_km(1,:), xyzNear_km(2,:), ...
            xyzNear_km(3,:));
        [rFar_km, thetaFar, phiFar] = xyz2sph(xyzFar_km(1,:), xyzFar_km(2,:), xyzFar_km(3,:));
        
        % Evaluate field differences in default model
        [MagModel, CsheetModel, ~, ~, ~] = GetModelOpts(planet, 0, 0);
        [BvecNear, ~, ~] = MagFldParent(planet, rNear_km, thetaNear, phiNear, MagModel, ...
            CsheetModel);
        [BvecFar, ~, ~] = MagFldParent(planet, rFar_km, thetaFar, phiFar, MagModel, CsheetModel);
        Bdiff = BvecNear - BvecFar;
        BdiffMag = sqrt(Bdiff(1,:).^2 + Bdiff(2,:).^2 + Bdiff(3,:).^2);
        maxBcompDiff = max(BdiffMag, [], 2);
        disp([moon ':'])
        disp(maxBcompDiff)
    end

end