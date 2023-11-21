function rhoMP_km = GetMPsurfSM1996(Rss_km, xyzPSMdipO_km)
% Determine paraboloid magnetopause surface shape based on Schulz and McNab (1996)
%
% Determines shape of magnetopause screening surface current shape from the model detailed by Schulz
% and McNabb (1996) https://doi.org/10.1029/95JA02987. This is a tanh paraboloid with an input
% sub-solar magnetopause standoff distance and dipole offset.
%
% Parameters
% ----------
% Rss_km : double, 1xN
%   Sub-solar point magnetopause standoff distance in km.
% xyzPSMdipO_km : double, 3xN
%   Dipole offset vector in planet--solar--magnetospheric coordinates in km.
%
% Returns
% -------
% rhoMP_km : double, 1xN
%   Distance from planet--Sun line in km, i.e. :math:`\rho` in dipole--solar--zenith cylindrical
%   coordinates.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rhoMP_km = 1.7696.*Rss_km .* tanh(0.61338./Rss_km .* (Rss_km - xyzPSMdipO_km(1,:))).^(1/2.280);
end
