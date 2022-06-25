function [rMinMoon_km, rJup_km] = GetMoonDist(sc, parentName, ets)

    switch parentName
        case 'Jupiter'
            moons = ["Io", "Europa", "Ganymede", "Callisto"];
        case 'Saturn'
            moons = ["Enceladus", "Rhea", "Dione", "Mimas", "Titan"];
        case 'Uranus'
            moons = ["Miranda", "Ariel", "Umbriel", "Titania", "Oberon"];
        case 'Neptune'
            moons = ["Triton"];
    end
    nMoons = length(moons);
    spkParent = upper(parentName);
    spkS3 = ['IAU_' spkParent];
    [scPos_km, ~] = cspice_spkpos(upper(sc), ets, spkS3, 'NONE', spkParent);
    spkMoonPos_km = zeros(nMoons, 3, length(ets));
    rMoons = zeros(nMoons, length(ets));
    for i=1:nMoons
        spkMoon = upper(char(moons(i)));
        [spkMoonPos_km(i,:,:), ~] = cspice_spkpos(spkMoon, ets, spkS3, 'NONE', spkParent);
        xDiff = squeeze(spkMoonPos_km(i,1,:))' - scPos_km(1,:);
        yDiff = squeeze(spkMoonPos_km(i,2,:))' - scPos_km(2,:);
        zDiff = squeeze(spkMoonPos_km(i,3,:))' - scPos_km(3,:);
        rMoons(i,:) = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2);
    end
    rMinMoon_km = squeeze(min(rMoons, [], 1));
    rJup_km = sqrt(scPos_km(1,:).^2 + scPos_km(2,:).^2 + scPos_km(3,:).^2);
end