function CalcBexc
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