function xMP_km = GetMPsurfAB2005(Rss_km, xyzPSM_km)
% Determine paraboloid magnetopause surface shape based on Alexeev and Belenkaya (2005).
%
% Determines shape of magnetopause screening surface current shape from the model detailed by
% Alexeev and Belenkaya (2005) https://doi.org/10.5194/angeo-23-809-2005. This is a simple
% paraboloid with an input sub-solar magnetopause standoff distance.
%
% Parameters
% ----------
% Rss_km : double, 1xN
%   Sub-solar point magnetopause standoff distance in km.
% xyzPSM_km : double, 3xN
%   Cartesian coordinates of measurement locations in planet--solar--magnetospheric frame in km.
%
% Returns
% -------
% xMP_km : double, 1xN
%   Distance from body planet--solar--magnetospheric yz-plane to magnetopause surface for each (y,z)
%   point, i.e. the x coordinate of the magnetopause in the PSM frame.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xMP_km = Rss_km - (xyzPSM_km(2,:).^2 + xyzPSM_km(3,:).^2)/2 ./ Rss_km;
end
