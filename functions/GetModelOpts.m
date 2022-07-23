function [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt)
    switch parentName
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
            
            AB2005 = 'AB2005'; % Alexeev and Belenkaya (2005) magnetopause model
            E1992a90 = 'Engle1992alpha90'; % Engle (1992) no-tilt magnetopause field model
            E1992a0 = 'Engle1992alpha0'; % Engle (1992) sunward-tilt magnetopause field model
            E1992a180 = 'Engle1992alpha180'; % Engle (1992) anti-sunward-tilt magnetopause field model
            B1994 = 'Bode1994'; % Engle (1992) magnetopause field model with Bode (1994) time-dependent coefficients
            
            MPmodel = B1994;

            if opt == 0; opt = 7; end  % Set default to JRM33 + C2020
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
                    magModelDescrip = 'Khurana & Schwarzl 2005';
                    fEnd = 'KS2005';
                case 4
                    MagModel = JRM09;
                    CsheetModel = C2020sheet;
                    magModelDescrip = 'JRM09 + C2020';
                    fEnd = 'JRM09C2020';
                case 5
                    MagModel = JRM09;
                    CsheetModel = C1981sheet;
                    magModelDescrip = 'Vance 2021 (JRM09 + C1981)';
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
                    disp(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            
        case 'Saturn'
            % opt 1
            B2010 = 'MagFldSaturnBurton2010';
            B2010sheet = 'SphericalHarmonic';
            % opt 2
            Cassini11 = 'MagFldSaturnCassini11';
            Cassini11sheet = 'Cassini11';
            
            MPmodel = 'None'; % No magnetopause model implemented

            if opt == 0; opt = 2; end  % Set default to Cassini 11
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
                    disp(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            
        case 'Uranus'
            % opt 1
            AH5 = 'MagFldUranusAH5';
            AH5sheet = 'None';
            MPmodel = 'None'; % No magnetopause model implemented

            if opt == 0; opt = 1; end  % Set default to AH5
            switch opt
                case 1
                    MagModel = AH5;
                    CsheetModel = AH5sheet;
                    magModelDescrip = 'AH5';
                    fEnd = 'AH5';
                otherwise
                    disp(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            
        case 'Neptune'
            % opt 1
            O8 = 'MagFldNeptuneO8';
            O8sheet = 'None';
            MPmodel = 'None'; % No magnetopause model implemented
            
            if opt == 0; opt = 1; end  % Set default to O8
            switch opt
                case 1
                    MagModel = O8;
                    CsheetModel = O8sheet;
                    magModelDescrip = 'O8';
                    fEnd = 'O8';
                otherwise
                    disp(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
    end
end