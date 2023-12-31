function CalcBexc
% Calculate excitation moments for all moons using selected eras and default models, in addition to
% particular selections for Europa and Callisto to allow for comparison to Seufert et al. (2011).

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PlanetMag('Europa', 'Galileo', 'SPRH', 1, 0, 0, 0, 0, 1e6, 6); % Seufert comparison
    PlanetMag('Callisto', 'Galileo', 'SPRH', 1, 0, 0, 0, 0, 1e6, 6); % Seufert comparison
    
    % Excitation moments for each moon for results
    PlanetMag('Io', 'Juno', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Europa', 'Juno', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Ganymede', 'Juno', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Callisto', 'Galileo', 'IAU', 1, 0, 0, 0, 0, 1e6, 6);

    PlanetMag('Mimas', 'Cassini', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Enceladus', 'Cassini', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Dione', 'Cassini', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Rhea', 'Cassini', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Titan', 'Cassini', 'IAU', 1, 0, 0, 0, 0, 1e6);

    PlanetMag('Miranda', 'Voyager', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Ariel', 'Voyager', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Umbriel', 'Voyager', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Titania', 'Voyager', 'IAU', 1, 0, 0, 0, 0, 1e6);
    PlanetMag('Oberon', 'Voyager', 'IAU', 1, 0, 0, 0, 0, 1e6);

    PlanetMag('Triton', 'Voyager', 'IAU', 1, 0, 0, 0, 0, 1e6);
end