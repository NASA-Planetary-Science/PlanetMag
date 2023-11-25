function ctimes = KS_ctimer(ets)
% Converts ephemeris times from TDB seconds relative to J2000 to ctime.
%
% MJS note: I don't understand why the reference is against 1/1/1966. The JS3(1965) epoch is
% 1965-01-01 at midnight (see https://doi.org/10.1029/GL004i002p00065), but the ctime columns in
% the trajectory files provided by K. Khurana for validation are consistent with a 1966 reference
% time.
%
% Parameters
% ----------
% ets : double, 1xN
%   Ephemeris times in TDB seconds relative to J2000.
%
% Returns
% -------
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J1965 = cspice_str2et('1966-01-01T00:00:00.000');
    adj = cspice_str2et('2000-01-01T12:00:00.000') + 37.5;
    ctimes = ets - J1965 + adj;
end
