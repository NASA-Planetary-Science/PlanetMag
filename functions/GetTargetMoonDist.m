function r_RM = GetTargetMoonDist(sc, moonName, parentName, ets)

    [~, ~, RM_km, ~, ~, ~, ~, ~, ~, ~] = GetBodyParams(moonName);
    spkParent = upper(parentName);
    spkS3 = ['IAU_' spkParent];
    [scPos_km, ~] = cspice_spkpos(upper(sc), ets, spkS3, 'NONE', spkParent);
    spkMoon = upper(moonName);
    [spkMoonPos_km, ~] = cspice_spkpos(spkMoon, ets, spkS3, 'NONE', spkParent);
    xDiff = spkMoonPos_km(1,:) - scPos_km(1,:);
    yDiff = spkMoonPos_km(2,:) - scPos_km(2,:);
    zDiff = spkMoonPos_km(3,:) - scPos_km(3,:);
    r_RM = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2) / RM_km;
end