function [eps1, eps2, ByConst, BzConst] = KS_coeffsBIMF
% Coefficients describing the interplanetary magnetic field (IMF) in calculations for the Khurana and
% Schwarzl (2005) model.
%
% Returns
% -------
% eps1, eps2 : double
%   Scaling factors for IMF parameters.
% ByConst, BzConst : double
%   Magnetic field components to use for the IMF.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eps1 = 0.068;
	eps2 = 0.554;
    ByConst = 0.0;
	BzConst = 0.001;
end
