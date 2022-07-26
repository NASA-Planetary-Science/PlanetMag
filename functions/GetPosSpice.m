function [r_km, theta_rad, phi_rad, xyz_km, S3coords] = GetPosSpice(moonName, parentName, t_h)
    t = t_h * 3600;
    spkParent = upper(parentName);
    spkMoon = upper(moonName);
    if strcmp(parentName, 'Uranus')
        S3coords = 'US3';
    else
        S3coords = ['IAU_' spkParent];
    end
    [xyz_km, ~] = cspice_spkpos(spkMoon, t, S3coords, 'NONE', spkParent);
    r_km = sqrt(xyz_km(1,:).^2 + xyz_km(2,:).^2 + xyz_km(3,:).^2);
    theta_rad = acos(xyz_km(3,:) ./ r_km);
    phi_rad = atan2(xyz_km(2,:), xyz_km(1,:));
end