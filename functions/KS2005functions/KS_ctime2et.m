function ets = KS_ctime2et(ctimes)
% Converts times from ctime to TDB seconds relative to J2000.
%
% Parameters
% ----------
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
%
% Returns
% -------
% ets : double, 1xN
%   Ephemeris times in TDB seconds relative to J2000.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J1965 = cspice_str2et('1966-01-01T00:00:00.000');
    adj = cspice_str2et('2000-01-01T12:00:00.000') + 37.5;
    ets = ctimes + J1965 - adj;
end
