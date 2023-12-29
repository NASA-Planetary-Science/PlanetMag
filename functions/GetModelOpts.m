function [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(planet, opt, MPopt)
% Convert from an index number to a compatible combination of planetary field models.
%
% For each planet, only specific internal and external magnetic field models are mutually
% compatible, as some are derived using parallel assumptions or outputs from other models. For
% ease of implementation, this function provides an index-based lookup to select from among valid
% model combinations.
%
% Parameters
% ----------
% planet : char, 1xC
%   Name of planet for which to return model combination parameters.
% opt : int
%   Index number for magnetospheric field model. 0 selects the default model, typically the most
%   recently, or most widely accepted, published model combination. Currently implemented model
%   combinations are:
%
%       - ``Earth``:
%           - ``1`` (or ``0``): IGRF-13 from Alken et al. (2021)
%             https://doi.org/10.1186/s40623-020-01288-x
%       - ``Jupiter``:
%           - ``1``: VIP4 model from Connerney et al. (1998) https://doi.org/10.1029/97JA03726 +
%             Connerney et al. (1981) current sheet model https://doi.org/10.1029/JA086iA10p08370
%           - ``2``: O6 model from Connerney (1992) https://core.ac.uk/download/pdf/83644007.pdf +
%             Khurana (1997) current sheet model https://doi.org/10.1029/97JA00563
%           - ``3``: VIP4 model with O6 model dipole moment orientation + associated current sheet
%             model, together from Khurana and Schwarzl (2005) https://doi.org/10.1029/2004JA010757
%           - ``4``: JRM09 model from Connerney et al. (2018) https://doi.org/10.1002/2018GL077312
%             + Connerney et al. (2020) current sheet model https://doi.org/10.1029/2020JA028138
%           - ``5``: JRM09 model from Connerney et al. (2018) https://doi.org/10.1002/2018GL077312
%             + Connerney et al. (1981) current sheet model https://doi.org/10.1029/JA086iA10p08370
%             This model combination is as applied in Vance et al. (2021)
%             https://doi.org/10.1029/2020JE006418
%           - ``6``: VIP4 model from Connerney et al. (1998) https://doi.org/10.1029/97JA03726 +
%             Khurana (1997) current sheet model https://doi.org/10.1029/97JA00563 This model
%             combination is as applied in Seufert et al. (2011)
%             https://doi.org/10.1016/j.icarus.2011.03.017
%           - ``7`` (or ``0``): JRM33 model from Connerney et al. (2021)
%             https://doi.org/10.1029/2021JE007138 + Connerney et al. (2020) current sheet model
%             https://doi.org/10.1029/2020JA028138.
%       - ``Saturn``:
%           - ``1``: Burton et al. (2010) model with degree-1 external field
%             https://doi.org/10.1029/2010GL045148
%           - ``2`` (or ``0``): Cassini 11 + corresponding current sheet model
%             https://doi.org/10.1126/science.aav6732
%       - ``Uranus``:
%           - ``1``: Q3 model with degree-1 external field from Connerney et al. (1987)
%             https://doi.org/10.1029/JA092iA13p15329
%           - ``2`` (or ``0``): AH5 model from Herbert (2009) https://doi.org/10.1029/2009JA014394
%       - ``Neptune``:
%           - ``1`` (or ``0``): O8 model from Connerney et al. (1992)
%             https://doi.org/10.1016/0273-1177(92)90394-D
%
% MPopt : int
%   Index number for magnetopause current magnetic field model. 0 selects the default option. If
%   any value is passed that does not correspond to a valid model index, no magnetopause current
%   magnetic field will be modeled. Currently implemented model options are:
%
%       - ``Earth``: No magnetopause current magnetic field model is implemented.
%       - ``Jupiter``:
%           - ``1``: Alexeev and Belenkaya (2005) model https://doi.org/10.5194/angeo-23-809-2005
%           - ``2`` (or ``0``): Engle (1992) model with no sunward or anti-sunward tilt
%             https://doi.org/10.1029/92JA02036
%           - ``3``: Engle (1992) model with sunward tilt
%             https://doi.org/10.1029/92JA02036
%           - ``4``: Engle (1992) model with anti-sunward tilt
%             https://doi.org/10.1029/92JA02036
%           - ``5``: Engle (1992) model  with Bode (1994) time-dependent
%             coefficients https://apps.dtic.mil/sti/citations/ADA284857
%       - ``Saturn``: No magnetopause current magnetic field model is implemented.
%       - ``Uranus``:
%           - ``1``: Q3mp model from Herbert (2009) https://doi.org/10.1029/2009JA014394
%           - ``2`` (or ``0``): Arridge and Eggington (2021) box harmonic model
%             https://doi.org/10.1016/j.icarus.2021.114562
%       - ``Neptune``:
%           - ``1`` (or ``0``): Schulz and McNabb (1996) model https://doi.org/10.1029/95JA02987
%
% Returns
% -------
% MagModel : char, 1xD
%   Planet intrinsic field model.
% CsheetModel : char, 1xE
%   Current sheet magnetic field model to apply, if applicable.
% MPmodel : char, 1xF
%   Magnetopause current magnetic field model to apply, if applicable.
% magModelDescrip : char, 1xG
%   Text description of magnetic field model combination for plot labels.
% fEnd : char, 1xH
%   File name descriptor to append for model-specific outputs, such as excitation moments.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch planet
        case 'Earth'
            % opt 1
            IGRF13 = 'MagFldEarthIGRF13';
            IGRF13sheet = 'None';
            
            if opt == 0; opt = 1; end % Set default to IGRF13
            switch opt
                case 1
                    MagModel = IGRF13;
                    CsheetModel = IGRF13sheet;
                    magModelDescrip = 'IGRF13';
                    fEnd = 'IGRF13';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ...
                        ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end

            MPmodel = 'None'; % No magnetopause model implemented
            MPend = 'noMP';

        case 'Jupiter'
            % opt 1
            VIP4 = 'MagFldJupiterVIP4';
            C1981sheet = 'Connerney1981';
            % opt 2
            O6 = 'MagFldJupiterGSFCO6';
            K1997sheet = 'Khurana1997';
            % opt 3
            KhuranaJup = 'VIP4 with O6 orientation';
            K2005sheet = 'KS2005';
            % opt 4
            JRM09 = 'MagFldJupiterJRM09';
            C2020sheet = 'Connerney2020';
            % opt 5 is Vance et al. 2021 combo, JRM09 + C1981
            % opt 6 is Seufert et al. 2011 combo, VIP4 + K1997
            % opt 7 (includes C2020 current sheet)
            JRM33 = 'MagFldJupiterJRM33';
            
            % MPopt 1: Alexeev and Belenkaya (2005) magnetopause model
            AB2005 = 'AB2005'; 
            % MPopt 2: Engle (1992) no-tilt magnetopause field model
            E1992a90 = 'Engle1992alpha90';
            % MPopt 3: Engle (1992) sunward-tilt magnetopause field model
            E1992a0 = 'Engle1992alpha0'; 
            % MPopt 4: Engle (1992) anti-sunward-tilt magnetopause field model
            E1992a180 = 'Engle1992alpha180';
            % MPopt 5: Engle (1992) magnetopause field model with Bode (1994) time-dependent
            % coefficients
            B1994 = 'Bode1994';

            if opt == 0; opt = 7; end % Set default to JRM33 + C2020
            switch opt
                case 1
                    MagModel = VIP4;
                    CsheetModel = C1981sheet;
                    magModelDescrip = 'VIP4 + C1981';
                    fEnd = 'VIP4C1981';
                case 2
                    MagModel = O6;
                    CsheetModel = K1997sheet;
                    magModelDescrip = 'O6 + K1997';
                    fEnd = 'O6K1997';
                case 3
                    MagModel = KhuranaJup;
                    CsheetModel = K2005sheet;
                    magModelDescrip = 'KS2005';
                    fEnd = 'KS2005';
                    MPopt = -1; % Prevent evalaution of an MP model because KS2005 already has one
                case 4
                    MagModel = JRM09;
                    CsheetModel = C2020sheet;
                    magModelDescrip = 'JRM09 + C2020';
                    fEnd = 'JRM09C2020';
                case 5
                    MagModel = JRM09;
                    CsheetModel = C1981sheet;
                    magModelDescrip = 'JRM09 + C1981';
                    fEnd = 'JRM09C1981';
                case 6
                    MagModel = VIP4;
                    CsheetModel = K1997sheet;
                    magModelDescrip = 'VIP4 + K1997';
                    fEnd = 'VIP4K1997';
                case 7
                    MagModel = JRM33;
                    CsheetModel = C2020sheet;
                    magModelDescrip = 'JRM33 + C2020';
                    fEnd = 'JRM33C2020';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ...
                        ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            if MPopt == 0; MPopt = 2; end % Set default to Engle (1992) no-tilt model
            switch MPopt
                case 1
                    MPmodel = AB2005;
                    MPend = 'AB2005';
                case 2
                    MPmodel = E1992a90;
                    MPend = 'E1992a90';
                case 3
                    MPmodel = E1992a0;
                    MPend = 'E1992a0';
                case 4
                    MPmodel = E1992a180;
                    MPend = 'E1992a180';
                case 5
                    MPmodel = B1994;
                    MPend = 'B1994';
                otherwise
                    MPmodel = 'None';
                    MPend = 'noMP';
            end
            
            
        case 'Saturn'
            % opt 1
            B2010 = 'MagFldSaturnBurton2010';
            B2010sheet = 'SphericalHarmonic';
            % opt 2
            Cassini11 = 'MagFldSaturnCassini11';
            Cassini11sheet = 'Cassini11';

            if opt == 0; opt = 2; end % Set default to Cassini 11
            switch opt
                case 1
                    MagModel = B2010;
                    CsheetModel = B2010sheet;
                    magModelDescrip = 'Burton et al. 2010';
                    fEnd = 'B2010';
                case 2
                    MagModel = Cassini11;
                    CsheetModel = Cassini11sheet;
                    magModelDescrip = 'Cassini 11 field + sheet';
                    fEnd = 'Cassini11';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ...
                        ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            MPmodel = 'None'; % No magnetopause model implemented
            MPend = 'noMP';
            
            
        case 'Uranus'
            % opt 1
            Q3 = 'MagFldUranusQ3';
            Q3sheet = 'SphericalHarmonic';
            % opt 2
            AH5 = 'MagFldUranusAH5';
            AH5sheet = 'None';
            
            % MPopt 1: Dipole-mirror magnetopause model used in Herbert
            % (2009) to derive AH5 model (based on SM1996)
            Q3mp = 'Q3mp';
            % MPopt 2: Box harmonic magnetopause model from Arridge and
            % Eggington (2021)
            AE2021 = 'AE2021';

            if opt == 0; opt = 2; end % Set default to AH5
            switch opt
                case 1
                    MagModel = Q3;
                    CsheetModel = 'None'; % Only use uniform external field if noMP
                    magModelDescrip = 'Q3';
                    fEnd = 'Q3';
                case 2
                    MagModel = AH5;
                    CsheetModel = AH5sheet;
                    magModelDescrip = 'AH5';
                    fEnd = 'AH5';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ...
                        ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            if MPopt == 0; MPopt = 2; end % Set default to AE2021 model
            switch MPopt
                case 1
                    MPmodel = Q3mp;
                    MPend = 'Q3mp';
                case 2
                    MPmodel = AE2021;
                    MPend = 'AE2021';
                otherwise
                    MPmodel = 'None';
                    MPend = 'noMP';
                    if opt == 1
                        CsheetModel = Q3sheet;
                    end
            end
            
            
        case 'Neptune'
            % opt 1
            O8 = 'MagFldNeptuneO8';
            O8sheet = 'None';
            
            % MPopt 1: Dipole-mirror magnetopause model described
            % in Schulz and McNab (1996)
            SM1996 = 'SM1996';
            
            if opt == 0; opt = 1; end % Set default to O8
            switch opt
                case 1
                    MagModel = O8;
                    CsheetModel = O8sheet;
                    magModelDescrip = 'O8';
                    fEnd = 'O8';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ...
                        ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            if MPopt == 0; MPopt = 1; end % Set default to SM1996 model
            switch MPopt
                case 1
                    MPmodel = SM1996;
                    MPend = 'SM1996';
                otherwise
                    MPmodel = 'None';
                    MPend = 'noMP';
            end
            
    end
    
    if ~strcmp(magModelDescrip, 'None')
        magModelDescrip = [magModelDescrip ' + ' MPend];
        fEnd = [fEnd MPend];
    end
end
