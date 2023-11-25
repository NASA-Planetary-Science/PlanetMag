function etime = KS_etimer(ctimes)
% Converts times from ctime to etime, which accounts for leap seconds.
%
% Parameters
% ----------
% ctimes : double, 1xN
%   Seconds past midnight Jan 1, 1966. See ctimer for more details.
%
% Returns
% -------
% etime : double, 1xN
%   ctimes adjusted for leap seconds.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tcor = zeros(size(ctimes));
    Tcor(ctimes >= 189302400.000) = Tcor(ctimes >= 189302400.000) + 10;
    Tcor(ctimes >= 205027200.000) = Tcor(ctimes >= 205027200.000) + 1;
    Tcor(ctimes >= 220924800.000) = Tcor(ctimes >= 220924800.000) + 1;
    Tcor(ctimes >= 252460800.000) = Tcor(ctimes >= 252460800.000) + 1;
    Tcor(ctimes >= 283996800.000) = Tcor(ctimes >= 283996800.000) + 1;
    Tcor(ctimes >= 315532800.000) = Tcor(ctimes >= 315532800.000) + 1;
    Tcor(ctimes >= 347155200.000) = Tcor(ctimes >= 347155200.000) + 1;
    Tcor(ctimes >= 378691200.000) = Tcor(ctimes >= 378691200.000) + 1;
    Tcor(ctimes >= 410227200.000) = Tcor(ctimes >= 410227200.000) + 1;
    Tcor(ctimes >= 441763200.000) = Tcor(ctimes >= 441763200.000) + 1;
    Tcor(ctimes >= 489024000.000) = Tcor(ctimes >= 489024000.000) + 1;
    Tcor(ctimes >= 520560000.000) = Tcor(ctimes >= 520560000.000) + 1;
    Tcor(ctimes >= 552096000.000) = Tcor(ctimes >= 552096000.000) + 1;
    Tcor(ctimes >= 615254400.000) = Tcor(ctimes >= 615254400.000) + 1;
    Tcor(ctimes >= 694224000.000) = Tcor(ctimes >= 694224000.000) + 1;
    Tcor(ctimes >= 757382400.000) = Tcor(ctimes >= 757382400.000) + 1;
    Tcor(ctimes >= 788918400.000) = Tcor(ctimes >= 788918400.000) + 1;
    Tcor(ctimes >= 836179200.000) = Tcor(ctimes >= 836179200.000) + 1;
    Tcor(ctimes >= 867715200.000) = Tcor(ctimes >= 867715200.000) + 1;
    Tcor(ctimes >= 899251200.000) = Tcor(ctimes >= 899251200.000) + 1;
    Tcor(ctimes >= 946684800.000) = Tcor(ctimes >= 946684800.000) + 1;
    Tcor(ctimes >= 993945600.000) = Tcor(ctimes >= 993945600.000) + 1;
    Tcor(ctimes >= 1041379200.00) = Tcor(ctimes >= 1041379200.00) + 1;
    etime = ctimes + Tcor - 0.1072958367816D10;
end
