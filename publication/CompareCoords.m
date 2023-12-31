function CompareCoords
% Compare phi-Omega frames to IAU frames for moons for which these frames are implemented (the
% Galilean moons). Prints the angles between x and z axes.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    moons = [
        "Io"
        "Europa"
        "Ganymede"
        "Callisto"
        ];
    
    % Iterate over Galilean moons
    for i=1:length(moons)
        moon = char(moons(i));
        spkMoon = upper(['IAU_' moon]);
        spkmPhiO = upper([moon '_PHI_O']);

        LoadSpice(moon,'None');
        [~, ~, RmoonEq_km, ~, ~, ~, ~, Tmoon_s, ~, ~] = GetBodyParams(moon);

        % Evaluate every 10 seconds over 2 orbital periods
        ets = 0:10:(2*Tmoon_s);
        orbFrac = ets/Tmoon_s;
        
        % Transform from IAU frame to PhiO frame for x and z axes
        zees = zeros(size(ets));
        [xmPhiO_km, ymPhiO_km, zmPhiO_km] = RotateVecSpice(RmoonEq_km, 0, 0, ets, spkmPhiO, ...
            spkMoon);
        x_xyzIAU_km = [zees + RmoonEq_km; zees; zees];
        x_xyzmPhiO_km = [xmPhiO_km; ymPhiO_km; zmPhiO_km];
        xAng_deg = rad2deg(cspice_vsep(x_xyzmPhiO_km, x_xyzIAU_km)) - 90;
        [xmPhiO_km, ymPhiO_km, zmPhiO_km] = RotateVecSpice(0, 0, RmoonEq_km, ets, spkmPhiO, ...
            spkMoon);
        z_xyzIAU_km = [zees; zees; zees + RmoonEq_km];
        z_xyzmPhiO_km = [xmPhiO_km; ymPhiO_km; zmPhiO_km];
        zAng_deg = rad2deg(cspice_vsep(z_xyzmPhiO_km, z_xyzIAU_km));

        % Plot the differences
        windowName = ['x axis difference for ' moon];
        titleInfo = ['Angle between IAU and ' moon(1) 'PhiO x axes'];
        xInfo = 'Fraction of orbital period past J2000';
        yInfo = 'x axis separation (deg)';
        fName = ['IAU_' moon(1) 'PhiO_xDiff'];
        PlotGeneric(orbFrac, xAng_deg, [], windowName, titleInfo, xInfo, yInfo, fName);

        windowName = ['z axis difference for ' moon];
        titleInfo = ['Angle between IAU and ' moon(1) 'PhiO z axes'];
        yInfo = 'z axis separation (deg)';
        fName = ['IAU_' moon(1) 'PhiO_zDiff'];
        PlotGeneric(orbFrac, zAng_deg, [], windowName, titleInfo, xInfo, yInfo, fName);
    end

end