function parentName = LoadSpice(moonName, sc, spiceDir)
% Load SPICE kernels needed to determine positions for an input spacecraft and target moon.
%
% Loads into the SPICE kernel pool all those kernels needed to evaluate spacecraft and moon
% positions relative to the parent planet.
%
% Parameters
% ----------
% moonName : char, 1xC
%   Name of target body. Currently implemented options are:
%   
%       -``Moon``
%       -``Io``
%       -``Europa``
%       -``Ganymede``
%       -``Callisto``
%       -``Mimas``
%       -``Enceladus``
%       -``Dione``
%       -``Rhea``
%       -``Titan``
%       -``Miranda``
%       -``Ariel``
%       -``Umbriel``
%       -``Titania``
%       -``Oberon``
%       -``Triton``
%
%   and their parent planets. The same binary kernels are loaded for any body within a particular
%   planetary system, as these kernels contain information for all listed satellites and the
%   planet.
% sc : char, 1xD
%   Spacecraft name for which trajectories will be loaded from relevant binary SPICE kernels.
%   Currently implemented options are:
%
%       -``Swarm``
%       -``Galileo``
%       -``Juno``
%       -``Voyager 1``
%       -``Voyager 2``
%
% spiceDir : char, 1xE, default='spice'
%   Directory where SPICE kernel files are located.
%
% Returns
% -------
% parentName : char, 1xF
%   Name of the parent body.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('sc', 'var'); sc = ''; end
    if ~exist('spiceDir', 'var'); spiceDir = 'spice'; end
    leapseconds = 'naif0012.tls';
    frames = 'custom_frames_v01.tf';
    
    if any(strcmp(moonName, {'Moon', 'Earth'}))
        parentName = 'Earth';
        generic = 'de430.bsp';
    elseif any(strcmp(moonName, {'Io', 'Europa', 'Ganymede', 'Callisto', 'Jupiter'}))
        parentName = 'Jupiter';
        generic = 'jup365.bsp';
        %generic = 'jup310.bsp';
    elseif any(strcmp(moonName, {'Mimas', 'Enceladus', 'Dione', 'Rhea', 'Titan', 'Saturn'}))
        parentName = 'Saturn';
        generic = 'sat441.bsp';
    elseif any(strcmp(moonName, {'Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon', 'Uranus'}))
        parentName = 'Uranus';
        generic = 'ura111.bsp';
    elseif any(strcmp(moonName, {'Triton', 'Neptune'}))
        parentName = 'Neptune';
        generic = 'nep097.bsp';
    else
        error([moonName ' has no parent planet defined in LoadSpice.'])
    end
    
    switch(sc)
        case 'Swarm'
            planetData = 'pck00010.tpc';
            pckEarth = 'earth_latest_high_prec.bpc';
            spiceKernelList = { leapseconds planetData generic pckEarth frames };
        case 'Galileo'
            planetData = 'pck00010.tpc';
            galileoSpicePrime = fullfile(sc, 's980326a.bsp');
            galileoSpiceGEM = fullfile(sc, 's000131a.bsp');
            galileoSpiceGMM = fullfile(sc, 's030916a.bsp');
            spiceKernelList = { leapseconds planetData galileoSpicePrime galileoSpiceGEM galileoSpiceGMM generic frames };
        case 'Juno'
            planetData = 'pck00010.tpc';
            generic = 'jup380s.bsp';
            junoGanyFlyby = fullfile(sc, 'juno_rec_210513_210630_210707.bsp');
            spiceKernelList = { leapseconds planetData generic junoGanyFlyby frames };
        case 'Voyager 1'
            planetData = 'pck00010.tpc';
            switch(parentName)
                case 'Jupiter'
                    voyaFlyby = fullfile(sc, 'vgr1_jup230.bsp');
                case 'Saturn'
                    voyaFlyby = fullfile(sc, 'vgr1_sat337.bsp');
            end
            spiceKernelList = { leapseconds planetData generic voyaFlyby frames };
        case 'Voyager 2'
            planetData = 'pck00010.tpc';
            switch(parentName)
                case 'Jupiter'
                    voyaFlyby = fullfile(sc, 'vgr2_jup230.bsp');
                case 'Saturn'
                    voyaFlyby = fullfile(sc, 'vgr2_sat337.bsp');
                case 'Uranus'
                    voyaFlyby = fullfile(sc, 'vgr2.ura111.bsp');
                case 'Neptune'
                    voyaFlyby = fullfile(sc, 'vgr2_nep097.bsp');
            end
            spiceKernelList = { leapseconds planetData generic voyaFlyby frames };
            
        otherwise
            planetData = 'pck00010.tpc';
            spiceKernelList = { leapseconds planetData generic frames };
    end
    
    cspice_furnsh(fullfile(spiceDir, spiceKernelList));
end