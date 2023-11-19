function [g, h, G, H, PlanetEqRadius, Nmax, NmaxExt] = GetGaussCoeffs(planet, InternalFieldModel)
% Read spherical harmonic coefficients from formatted ASCII files on disk.
%
% Spherical harmonic components for each model are returned in gauss, as are the maximum degree
% implemented in the reported model and in P\lanetMag and the planet reference radius in km used in
% the spherical harmonic expansion, as reported in listed publication from which the coefficients
% have been collected.
%
% Parameters
% ----------
% planet : char, 1xC
%   Planet for which to load model coefficients. Options are the planets listed under
%   ``InternalFieldModel``, as well as the special name ``'PureHarmonic'``, which reads a register
%   file from disk that has a single g or h represented for a specific n,m pair, for testing and
%   validation purposes.
% InternalFieldModel : char, 1xD
%   Published spherical harmonic model for the desired planet. Currently implemented are:
%
%       - ``Earth``
%           - ``MagFldEarthIGRF13``: "International Geomagnetic Reference Field: the thirteenth
%             generation." https://doi.org/10.1186/s40623-020-01288-x
%       - ``Jupiter``
%           - ``MagFldJupiterVIP4``: "New models of Jupiter's magnetic field constrained by the Io
%             flux tube footprint." https://doi.org/10.1029/97JA03726
%             (Voyager--Io footprint--Pioneer degree 4 model.)
%           - ``MagFldJupiterGSFCO4``: "The main magnetic field of Jupiter."
%             https://doi.org/10.1029/JA081i016p02917 The P\lanetMag implementation is based on the
%             update to System III (1965) conventions by Connerney (1992):
%             https://core.ac.uk/download/pdf/83644007.pdf. The model is named GSFC O4 after the 
%             institution where it was developed (Goddard Spaceflight Center, GSFC) and that it
%             includes up to the octupole (degree 3) coefficients from a degree 4 expansion. 
%           - ``MagFldJupiterGSFCO6``: Coefficients as repoted in "Magnetic Fields of the Outer
%             Planets." https://doi.org/10.1007/s11214-009-9621-7, though originally described in
%             Connerney (1992): https://core.ac.uk/download/pdf/83644007.pdf. Also developed at
%             GSFC, this model includes the octupole coefficients from a degree 6 expansion, hence
%             the name O6.
%           - ``MagFldJupiterJRM09``: "A New Model of Jupiter's Magnetic Field From Juno's First
%             Nine Orbits." https://doi.org/10.1002/2018GL077312 The Jupiter Reference Model
%             after 9 orbits by Juno.
%           - ``MagFldJupiterJRM33``: "A New Model of Jupiter's Magnetic Field at the Completion of
%             Juno's Prime Mission." https://doi.org/10.1029/2021JE007138 The Jupiter Reference
%             Model after 33 orbits of Jupiter by Juno.
%       - ``Saturn``
%           - ``MagFldSaturnBurton2010``: "Saturn's internal planetary magnetic field."
%             https://doi.org/10.1029/2010GL045148
%           - ``MagFldSaturnCassini11``: "Saturn's magnetic field revealed by the Cassini Grand
%             Finale." https://doi.org/10.1126/science.aat5434, corrected values at
%             https://doi.org/10.1126/science.aav6732
%       - ``Uranus``
%           - ``MagFldUranusQ3``: "The magnetic field of Uranus."
%             https://doi.org/10.1029/JA092iA13p15329 The name Q3 comes from including only a
%             dipole and quadrupole moment in a degree-3 spherical harmonic expansion. This model
%             also includes degree-1 external field coefficients, i.e. a uniform background field.
%           - ``MagFldUranusAH5``: "Aurora and magnetic field of Uranus."
%             https://doi.org/10.1029/2009JA014394 AH5 stands for aurora hexadecapole L-shell = 5
%             model.
%       - ``Neptune``:
%           - ``MagFldNeptuneO8``: "The magnetic field of Neptune."
%             https://doi.org/10.1016/0273-1177(92)90394-D Named similarly to the Jupiter GSFC
%             models, the Neptune O8 model is a degree-8 expansion that retains only up to the
%             octupole moment.
%
% Returns
% -------
% g, h : double, (Nmax)x(Nmax+1)
%   Real spherical harmonic coefficients (Gauss coefficients) in Schmidt semi-normalized form for
%   the internal magnetic field of the planet in terms of surface field strength in gauss.
%   Coefficients are each referenced to a specific System III coordinate system: ITRF08 for Earth,
%   ULS for Uranus, NLS for Neptune, and IAU for Jupiter and Saturn. Refer to the listed
%   publications for more details.
% G, H : double, (NmaxExt)x(NmaxExt+1)
%   Real spherical harmonic coefficients (Gauss coefficients) in Schmidt semi-normalized form for
%   the external magnetic field applied to the planet in gauss.
% PlanetEqRadius : double
%   Planet reference (equatorial) radius used in reporting spherical harmonic coefficients from
%   each reference.
% Nmax : int
%   Maximum degree of spherical harmonic coefficients in the selected internal field model.
% NmaxExt : int
%   Maximum degree of spherical harmonic coefficients in the external field model associated with
%   the selected internal field model, if applicable. Set to zero if the model does not include
%   external field coefficients.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    coeffPath = './modelCoeffs/';
    nHeadLines = 2;
    G = 0; 
    H = 0;
    NmaxExt = 0;
    switch(planet)
        case 'Earth'
            switch(InternalFieldModel)
                case 'MagFldEarthIGRF13'
                    % IGRF13, International Geomagnetic Reference Field:
                    % the thirteenth generation. https://doi.org/10.1186/s40623-020-01288-x
                    % Interior field parameters (Schmidt semi-normalized coefficients) in Gauss
                    % referenced to ITRF08 coordinates, 1RE = 6,371.2 km

                    Nmax = 13; % order
                    PlanetEqRadius = 6371.2; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsEarthIGRF13g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsEarthIGRF13h.csv']), ',', nHeadLines, 0);

                otherwise
                    error(['Magnetic field model "' InternalFieldModel '" not recognized.'])
            end

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
                    % From Russell and Dougherty (2010) Jupiter GSFC O6 magnetic field model:
                    % https://doi.org/10.1007/s11214-009-9621-7
                    % Interior field parameters (Schmidt semi-normalized Coefficients) in Gauss
                    % referenced to JSIII (1965) coordinates, 1RJ = 71,372 km.
                    % Originally reported in Connerney (1992):
                    % https://core.ac.uk/download/pdf/83644007.pdf

                    Nmax = 3; % order
                    PlanetEqRadius = 71372; % km, as reported in the publication

                    g = dlmread(fullfile([coeffPath 'coeffsJupiterO6g.csv']), ',', nHeadLines, 0);
                    h = dlmread(fullfile([coeffPath 'coeffsJupiterO6h.csv']), ',', nHeadLines, 0);

                case 'MagFldJupiterJRM09'
                    % Connerney et al. (2018), A New Model of Jupiter's Magnetic Field From Juno's
                    % First Nine Orbits
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
                    % Dougherty et al. (2018),  Saturn's magnetic field revealed by the Cassini
                    % Grand Finale
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

        case 'PureHarmonic'
            % Evaluate a single g,h spherical harmonic as detailed in the appropriate files. For
            % testing and validation purposes.

            Nmax = 10;
            PlanetEqRadius = 1;
            [g, h] = deal(zeros(Nmax, Nmax+1));
            pureHarm = dlmread(fullfile([coeffPath 'coeffsPureHarmonic.csv']), ',', nHeadLines, 0);
            n = pureHarm(1);
            m = pureHarm(2);
            g(n,m+1) = pureHarm(3);
            h(n,m+1) = pureHarm(4);
            
        otherwise
            error(['Planet name "' planet '" not recognized. Options are Jupiter, Saturn, Uranus, Neptune.'])

    end
end