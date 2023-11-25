function [period_hr, omegaJ_deghr, X0, phip0, obliq, incl_solar, C, psi2, psi4] ...
    = KS_coeffsCsheetStruc(AS_CODED)
% Import coefficients describing the structure of the current sheet in the Khurana and Schwarzl
% (2005) model.
% 
% Imports coefficients used to determine the current sheet structure used in calculating current
% sheet location relative to System III equatorial plane.
% 
% Parameters
% ----------
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% period_hr : double
%   Planetary rotation period in hours.
% omegaJ_deghr : double
%   Angular velocity of planetary rotation in degrees per hour.
% X0 : double
%   Distance from JSO :math:`x=0` plane beyond which the current sheet is hinged due to solar wind
%   forcing.
% phip0 : double
%   Current sheet parameterization coefficient.
% obliq : double
%   Obliquity of Jupiter relative to the solar spin axis in radians.
% incl_solar : double
%   Orbital inclination of Jupiter relative to the solar spin equator in radians.
% C : double, 1x5
%   Current sheet parameterization coefficient.
% psi2, psi4 : double
%   Current sheet parameterization coefficients.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if AS_CODED
        period_hr = 9.927953;
        omegaJ_deghr = 36.26125;
    else
        % Currently accepted System III period is 9.92492 h, not 9.927953.
        % Use SPICE PCK file if loaded
        try
            PMvals = cspice_bodvcd(599, 'PM', 3);
            omegaJ_deghr = PMvals(2) / 24;
        catch
            omegaJ_deghr = 870.536 / 24;
        end
        period_hr = 360 / omegaJ_deghr;
    end
    X0 = -45.0;
    phip0 = 6.12611;
    % MJS: The following numbers, called "three" and "six" in K. Khurana's code, are the obliquity
    % and orbital inclination relative to the solar spin equator.
    obliq = 0.0558505360;
    incl_solar = 0.11170107;
    C = [0.005973, 5.114e-5, 1.59e-5, 0.313244, -0.366166];
    psi2 = -1.201201;
    psi4 = 2.2604522;
end
