function [r_km, latS3_deg, lonS3_deg, xyz_km] = GetPosSpice(moonName, parentName, t_h)
    t = t_h * 3600;
    spkParent = upper(parentName);
    spkMoon = upper(moonName);
    spkS3 = ['IAU_' spkParent];
    [xyz_km, ~] = cspice_spkpos(spkMoon, t, spkS3, 'NONE', spkParent);
    r_km = sqrt(xyz_km(1,:).^2 + xyz_km(2,:).^2 + xyz_km(3,:).^2);
    latS3_deg = asind(xyz_km(3,:) ./ r_km);
    lonS3_deg = mod(atan2d(xyz_km(2,:), xyz_km(1,:)), 360.0);
end