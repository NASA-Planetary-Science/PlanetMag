function parentName = LoadSpice(moonName, sc)
    if ~exist('sc', 'var'); sc = ''; end
    spiceDir = 'spice/';
    leapseconds = 'naif0012.tls';
    PhiO = 'clipper_dyn_v01_mod.tf';
    
    if any(strcmp(moonName, {'Io', 'Europa', 'Ganymede', 'Callisto'}))
        parentName = 'Jupiter';
        generic = 'jup365.bsp';
        %generic = 'jup310.bsp';
    elseif any(strcmp(moonName, {'Mimas', 'Enceladus', 'Dione', 'Rhea', 'Titan'}))
        parentName = 'Saturn';
        generic = 'sat427.bsp';
    elseif any(strcmp(moonName, {'Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon'}))
        parentName = 'Uranus';
        generic = 'ura111.bsp';
    elseif any(strcmp(moonName, {'Triton'}))
        parentName = 'Neptune';
        generic = 'nep095.bsp';
    else
        error([moonName ' has no parent planet defined in LoadSpice.'])
    end
    
    if strcmp(sc, 'Galileo')
        planetData = 'pck00006.tpc';
        galileoSpicePrime = 's980326a.bsp';
        galileoSpiceGEM = 's000131a.bsp';
        galileoSpiceGMM = 's030916a.bsp';
        spiceKernelList = { leapseconds planetData galileoSpicePrime galileoSpiceGEM galileoSpiceGMM generic PhiO };
    elseif strcmp(sc, 'Juno')
        planetData = 'pck00010.tpc';
        generic = 'jup380s.bsp';
        junoGanyFlyby = 'juno_rec_210513_210630_210707.bsp';
        spiceKernelList = { leapseconds planetData generic junoGanyFlyby PhiO };
    else
        planetData = 'pck00010.tpc';
        spiceKernelList = { leapseconds planetData generic PhiO };
    end
    
    cspice_furnsh(fullfile(strcat(spiceDir, spiceKernelList)));
end