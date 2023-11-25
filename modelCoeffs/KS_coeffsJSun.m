function [yrJup_s, omegaJ_degps, omegayr_radps, etime1, obliq, tan_ob, aa, bb, deltaPhi] ...
    = KS_coeffsJSun(AS_CODED)
% Import coefficients for Jupiter--Sun angle calculations.
% 
% Parameters
% ----------
% AS_CODED : bool, default=0
%   Whether to match the original Fortran code (true) or with increased precision and corrected
%   parameters (false). See MagFldJupiterKS2005 for more details.
%
% Returns
% -------
% yrJup_s : double
%   Jupiter's orbital period around the Sun in seconds.
% omegaJ_degps : double
%   Sidereal rotation rate of Jupiter in degrees per second.
% omegayr_radps : double
%   Angular speed of orbital motion for Jupiter in radians per second.
% etime1 : double
%   Reference time from which orbital/rotation parameters are referenced.
% obliq : double
%   Obliquity of Jupiter relative to the solar spin axis in radians.
% tan_ob : double
%   ``tan(obliq)``.
% aa, bb : double, 1x7
%   Solar phase angle parameterization coefficients.
% deltaPhi : double
%   Sun direction offset in degrees.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    aa = [0.14347029, ...
          3.1145815, ...
         -0.12025561, ...
          0.093909436, ...
         -0.39321884e-5, ...
          0.10194945e-3, ...
         -0.12799464];
    bb = [-4.5467523, ...
           3.1848875, ...
          -0.16329986, ...
          -0.09776818, ...
           0.17556527e-3, ...
          -0.01978317, ...
           44.55915];
    deltaPhi = 48.23012;

    if AS_CODED
        yrJup_s = 11.85652502 * 365.25 * 86400;
        omegaJ_degps = 870.536 / 86400;
        thepi = 3.1415927;
        obliq = 3.123 * thepi/180;
        tan_ob = 0.054560676;
    else
        % SPICE can be used to calculate this angle directly. Doing so would provide a much more
        % precise way of making these calculations consistently.
        yrJup_s = 4332.589 * 24 * 86400;
        % Use SPICE PCK if loaded
        try
            PMvals = cspice_bodvcd(599, 'PM', 3);
            omegaJ_degps = PMvals(2) / 86400;
        catch
            omegaJ_degps = 870.536 / 86400;
        end
        thepi = pi;
        obliq = 3.123 * thepi/180;
        tan_ob = tan(obliq);
    end
    
    omegayr_radps = 2*thepi / yrJup_s;
    etime1 = -8.25767955817479D8;
end
