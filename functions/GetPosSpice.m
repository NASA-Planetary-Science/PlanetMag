function [r_km, theta_rad, phi_rad, xyz_km, S3coords] = GetPosSpice(moonName, parentName, t_h, S3coords)
    t = t_h * 3600;
    spkParent = upper(parentName);
    spkMoon = upper(moonName);
    if ~exist('S3coords', 'var')
        switch(parentName)
            case 'Uranus'
                S3coords = 'ULS';
            case 'Neptune'
                S3coords = 'NLS';
            otherwise
                S3coords = ['IAU_' spkParent];
        end
    end
    [xyz_km, ~] = cspice_spkpos(spkMoon, t, S3coords, 'NONE', spkParent);
    r_km = sqrt(xyz_km(1,:).^2 + xyz_km(2,:).^2 + xyz_km(3,:).^2);
    theta_rad = acos(xyz_km(3,:) ./ r_km);
    phi_rad = atan2(xyz_km(2,:), xyz_km(1,:));
end