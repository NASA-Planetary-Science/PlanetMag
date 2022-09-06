function [g, h, G, H, PlanetEqRadius, Nmax, NmaxExt] = GetGaussCoeffs(planet, InternalFieldModel)

    coeffPath = './modelCoeffs/';
    nHeadLines = 2;
    G = 0; 
    H = 0;
    NmaxExt = 0;
    switch(planet)
        case 'Jupiter'
            switch(InternalFieldModel)
                case 'MagFldJupiterVIP4'
                    % Connerney1998, New models of Jupiter's magnetic field constrained by the Io flux tube footprint
                    % Interior field parameters (Schmidt semi-normalized coefficients) in Gauss
                    % referenced to JSIII (1965) coordinates, 1RJ = 71,323 km

                    Nmax = 4; % order
                    PlanetEqRadius = 71323; % km, as reported in the publication (table of coefficients).
                    % Note that the text of Connerney et al. (1998) includes a
                    % conflicting definition of RJ = 71,398 km.

                    g = dlmread(fullfile([coeffPath 'coeffsJupiterVIP4g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsJupiterVIP4h.csv']), ',', nHeadLines, 0);

                case 'MagFldJupiterGSFCO4'
                    % From Acuna and Ness (1976) Jupiter GSFC O4 magnetic field model
                    % Interior field parameters (Schmidt semi-normalized Coefficients) in Gauss
                    % referenced to JSIII (1957) coordinates, 1RJ = 71,372 km. 
                    % Note that as the original publication used SIII (1957) coordinates,
                    % these values were updated to SIII (1965) coordinates by Acuna 
                    % et al. (1983), as reported in Connerney (1992) along with the 
                    % O6 model: https://core.ac.uk/download/pdf/83644007.pdf

                    Nmax = 3; % order
                    PlanetEqRadius = 71372; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsJupiterO4g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsJupiterO4h.csv']), ',', nHeadLines, 0);

                case 'MagFldJupiterGSFCO6'
                    % From Russell and Dougherty (2010) Jupiter GSFC O6 magnetic field model
                    % Interior field parameters (Schmidt semi-normalized Coefficients) in Gauss
                    % referenced to JSIII (1965) coordinates, 1RJ = 71,372 km.
                    % Originally reported in Connerney (1992): https://core.ac.uk/download/pdf/83644007.pdf

                    Nmax = 3; % order
                    PlanetEqRadius = 71372; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsJupiterO6g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsJupiterO6h.csv']), ',', nHeadLines, 0);

                case 'MagFldJupiterJRM09'
                    % Connerney et al. (2018), A New Model of Jupiter's Magnetic Field From Juno's First Nine Orbits
                    % JRM09: Juno Reference Model through perijove 9
                    % 10 (Perijove) PJ1-9, JUPITER, Juno, JRM09, Rc=0.85, I20 md, 264ev,
                    % r<7Rj 12/19/2017. lists the degree/order of the expansion (10) 
                    % and identifying information

                    Nmax = 10; % order
                    PlanetEqRadius = 71492; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsJupiterJRM09g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsJupiterJRM09h.csv']), ',', nHeadLines, 0);

                case 'MagFldJupiterJRM33'
                    % Connerney2020, A New Model of Jupiter's Magnetic Field at the Completion of Juno's Prime Mission
                    % JRM09: Juno Reference Model through perijove 33, on 2021-04-15
                    % "The degree 1 coefficients of the JRM33 model describe a dipole with moment M = 4.177 G,
                    %  offset from the rotation axis by θd = 10.25° towards System III longitude of ϕd = 196.38°.
                    %  Differences with the previous JRM09 model dipole (ΔM = 0.007 G, Δθd = 0.06°, Δϕd = 0.23°)
                    %  are slight."
                    % 30  1  JUPITER RJ= 71,492. rcore=0.8, r=[2.,2.5], I30 E1, 620ev, Jun 28, 2021

                    Nmax = 10; % order
                    PlanetEqRadius = 71492; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsJupiterJRM33g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsJupiterJRM33h.csv']), ',', nHeadLines, 0);

                otherwise
                    error(['Magnetic field model "' InternalFieldModel '" not recognized.'])
            end

        case 'Saturn'
            switch(InternalFieldModel)
                case 'MagFldSaturnBurton2010'
                    % Burton et al. (2010), Saturn's internal planetary magnetic field
                    % Interior field parameters (Schmidt semi-normalized coefficients) in Gauss
                    % referenced to SIII coordinates, 1RS = 60,268 km

                    Nmax = 3; % order
                    PlanetEqRadius = 60268; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsSaturnBurton2010g.csv']), ',', nHeadLines, 0);
                    h = zeros(Nmax,Nmax+1);

                    NmaxExt = 1;
                    G = dlmread(fullfile([coeffPath 'coeffsSaturnBurton2010Gext.csv']), ',', nHeadLines, 0);
                    H = dlmread(fullfile([coeffPath 'coeffsSaturnBurton2010Hext.csv']), ',', nHeadLines, 0);

                case 'MagFldSaturnCassini11'
                    % Dougherty et al. (2018),  Saturn's magnetic field revealed by the Cassini Grand Finale
                    % Interior field parameters (Schmidt semi-normalized coefficients) in Gauss
                    % referenced to SIII coordinates, 1RS = 60,268 km

                    Nmax = 12; % maximum resolved order is 11
                    PlanetEqRadius = 60268; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsSaturnCassini11g.csv']), ',', nHeadLines, 0);
                    h = zeros(Nmax,Nmax+1);

                otherwise
                        error(['Magnetic field model "' InternalFieldModel '" not recognized.'])
            end

        case 'Uranus'
            switch(InternalFieldModel)
                case 'MagFldUranusQ3'
                    % Connerney et al. (1987), The magnetic field of Uranus
                    % Q3: Quadrupole model with some unresolved coefficients reported
                    % up to n = 3 for Uranus, from Voyager 2 data.

                    Nmax = 3; % maximum resolved order is 2
                    PlanetEqRadius = 25600; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsUranusQ3g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsUranusQ3h.csv']), ',', nHeadLines, 0);

                    NmaxExt = 1;
                    G = dlmread(fullfile([coeffPath 'coeffsUranusQ3Gext.csv']), ',', nHeadLines, 0);
                    H = dlmread(fullfile([coeffPath 'coeffsUranusQ3Hext.csv']), ',', nHeadLines, 0);

                case 'MagFldUranusAH5'
                    % Herbert (2009), Aurora and magnetic field of Uranus
                    % AH5: Aurora hexadecapole L=5 model with unresolved coefficients reported
                    % up to n = 4 for Uranus, from Voyager 2 data.

                    Nmax = 4; % maximum resolved order is 3
                    PlanetEqRadius = 25559; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsUranusAH5g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsUranusAH5h.csv']), ',', nHeadLines, 0);

                otherwise
                    error(['Magnetic field model "' InternalFieldModel '" not recognized.'])

            end

        case 'Neptune'
            switch(InternalFieldModel)
                case 'MagFldNeptuneO8'
                    % Connerney et al. (1991), The magnetic field of Neptune
                    % O8: Octupole model with (mostly unresolved) coefficients reported
                    % up to n = 8 for Neptune, from Voyager 2 data.

                    Nmax = 8; % maximum resolved order is 3
                    PlanetEqRadius = 24765; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsNeptuneO8g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsNeptuneO8h.csv']), ',', nHeadLines, 0);

                otherwise
                    error(['Magnetic field model "' InternalFieldModel '" not recognized.'])

            end

        otherwise
            error(['Planet name "' planet '" not recognized. Options are Jupiter, Saturn, Uranus, Neptune.'])

    end
end